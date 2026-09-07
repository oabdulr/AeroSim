#include "wing3d.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define PI_F 3.14159265358979f
#define NSEC 64 /* samples per surface for the section outline */

static void camber_at(const Wing3 *w, float xc, float *yc, float *dyc)
{
    float m = w->camber, p = w->camber_pos;
    if (m <= 0.f || p <= 0.f || p >= 1.f) { *yc = 0.f; *dyc = 0.f; return; }
    if (xc < p) {
        *yc  = m / (p * p) * (2.f * p * xc - xc * xc);
        *dyc = 2.f * m / (p * p) * (p - xc);
    } else {
        float q = (1.f - p) * (1.f - p);
        *yc  = m / q * ((1.f - 2.f * p) + 2.f * p * xc - xc * xc);
        *dyc = 2.f * m / q * (p - xc);
    }
}

static float thickness_at(const Wing3 *w, float xc)
{
    float t = w->thickness;
    return 5.f * t * (0.2969f * sqrtf(xc) - 0.1260f * xc - 0.3516f * xc * xc
                      + 0.2843f * xc * xc * xc - 0.1036f * xc * xc * xc * xc);
}

int wing3_section(const Wing3 *w, int flapped, float *px, float *py)
{
    float sxu[NSEC], syu[NSEC], sxl[NSEC], syl[NSEC];

    float d  = flapped ? w->flap_deg * PI_F / 180.f : 0.f;
    float xh = w->flap_hinge;
    float ych, dych;
    camber_at(w, xh, &ych, &dych);
    float cd = cosf(d), sd = sinf(d);

    for (int k = 0; k < NSEC; ++k) {
        float sf = (float)k / (float)(NSEC - 1);
        float xc = 0.5f * (1.f - cosf(PI_F * sf));
        float yc, dyc;
        camber_at(w, xc, &yc, &dyc);
        float yt = thickness_at(w, xc);
        float th = atanf(dyc);
        float sn = sinf(th), cs = cosf(th);
        float xu = xc - yt * sn, yu = yc + yt * cs;
        float xl = xc + yt * sn, yl = yc - yt * cs;
        if (d != 0.f && xu > xh) {
            float dx = xu - xh, dy = yu - ych;
            xu = xh + dx * cd + dy * sd;
            yu = ych - dx * sd + dy * cd;
        }
        if (d != 0.f && xl > xh) {
            float dx = xl - xh, dy = yl - ych;
            xl = xh + dx * cd + dy * sd;
            yl = ych - dx * sd + dy * cd;
        }
        sxu[k] = xu; syu[k] = yu;
        sxl[k] = xl; syl[k] = yl;
    }

    int n = 0;
    for (int k = NSEC - 1; k >= 0; --k) { px[n] = sxu[k]; py[n] = syu[k]; ++n; }
    for (int k = 1; k < NSEC; ++k)      { px[n] = sxl[k]; py[n] = syl[k]; ++n; }

    /* chord scale, AoA rotation about quarter chord, translate to grid */
    float a  = w->aoa_deg * PI_F / 180.f;
    float ca = cosf(a), sa = sinf(a);
    for (int k = 0; k < n; ++k) {
        float X = px[k] * w->chord - 0.25f * w->chord;
        float Y = py[k] * w->chord;
        px[k] = w->cx + X * ca + Y * sa;
        py[k] = w->cy - X * sa + Y * ca;
    }
    return n;
}

/* Nonzero-winding scanline fill of one section into one z-slice. */
static void fill_slice(Sim3 *s, uint8_t *mask, int z,
                       const float *px, const float *py, int n)
{
    float ymin = 1e9f, ymax = -1e9f;
    for (int k = 0; k < n; ++k) {
        if (py[k] < ymin) ymin = py[k];
        if (py[k] > ymax) ymax = py[k];
    }
    (void)s;
    int y0 = (int)floorf(ymin) - 1; if (y0 < 1) y0 = 1;
    int y1 = (int)ceilf(ymax) + 1;  if (y1 > NY3 - 2) y1 = NY3 - 2;

    for (int y = y0; y <= y1; ++y) {
        float xs[32]; int dir[32]; int m = 0;
        float fy = (float)y;
        for (int k = 0; k < n; ++k) {
            float ax = px[k], ay = py[k];
            float bx = px[(k + 1) % n], by = py[(k + 1) % n];
            if ((ay > fy) != (by > fy)) {
                float t = (fy - ay) / (by - ay);
                if (m < 32) { xs[m] = ax + t * (bx - ax); dir[m] = by > ay ? 1 : -1; ++m; }
            }
        }
        for (int a = 1; a < m; ++a) {
            float vx = xs[a]; int vd = dir[a]; int b = a - 1;
            while (b >= 0 && xs[b] > vx) { xs[b+1] = xs[b]; dir[b+1] = dir[b]; --b; }
            xs[b+1] = vx; dir[b+1] = vd;
        }
        int wind = 0;
        for (int k = 0; k + 1 < m; ++k) {
            wind += dir[k];
            if (wind != 0) {
                int xa = (int)ceilf(xs[k]);
                int xb = (int)floorf(xs[k + 1]);
                if (xa < 1) xa = 1;
                if (xb > NX3 - 2) xb = NX3 - 2;
                for (int x = xa; x <= xb; ++x) mask[IDX3(x, y, z)] = 1;
            }
        }
    }
}

void wing3_rasterize(const Wing3 *w, Sim3 *s)
{
    static uint8_t *news = NULL;
    if (!news) news = (uint8_t *)malloc(NCELLS3);
    memset(news, 0, NCELLS3);

    float pxc[WG_MAX_PTS], pyc[WG_MAX_PTS]; /* clean section  */
    float pxf[WG_MAX_PTS], pyf[WG_MAX_PTS]; /* flapped section */
    int nc = wing3_section(w, 0, pxc, pyc);
    int nf = wing3_section(w, 1, pxf, pyf);

    int z0 = (int)floorf(w->cz - 0.5f * w->span); if (z0 < 1) z0 = 1;
    int z1 = (int)ceilf (w->cz + 0.5f * w->span); if (z1 > NZ3 - 2) z1 = NZ3 - 2;
    float halfflap = 0.5f * w->flap_frac * w->span;

    for (int z = z0; z <= z1; ++z) {
        int flapped = fabsf((float)z - w->cz) <= halfflap;
        if (flapped) fill_slice(s, news, z, pxf, pyf, nf);
        else         fill_slice(s, news, z, pxc, pyc, nc);
    }

    for (int i = 0; i < NCELLS3; ++i) {
        if (s->solid[i] && !news[i]) {
            sim3_refresh_cell(s, i);
        } else if (!s->solid[i] && news[i]) {
            s->rho[i] = 1.f;
            s->ux[i] = s->uy[i] = s->uz[i] = 0.f;
        }
    }
    memcpy(s->solid, news, NCELLS3);
    sim3_update_careful_mask(s);
}

/* render mesh */

static void emit_tri(float **o, int *nv,
                     const float a[3], const float b[3], const float c[3])
{
    float e1[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
    float e2[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
    float n[3]  = { e1[1]*e2[2]-e1[2]*e2[1],
                    e1[2]*e2[0]-e1[0]*e2[2],
                    e1[0]*e2[1]-e1[1]*e2[0] };
    float l = sqrtf(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (l > 1e-12f) { n[0] /= l; n[1] /= l; n[2] /= l; }
    const float *v[3] = { a, b, c };
    for (int k = 0; k < 3; ++k) {
        (*o)[0] = n[0]; (*o)[1] = n[1]; (*o)[2] = n[2];
        (*o)[3] = v[k][0]; (*o)[4] = v[k][1]; (*o)[5] = v[k][2];
        *o += 6; *nv += 1;
    }
}

int wing3_mesh(const Wing3 *w, float *out, int maxverts)
{
    enum { NST = 48 }; /* spanwise stations */
    static float rx[NST + 1][2 * NSEC], ry[NST + 1][2 * NSEC];
    int nring = 0;

    float pxc[WG_MAX_PTS], pyc[WG_MAX_PTS], pxf[WG_MAX_PTS], pyf[WG_MAX_PTS];
    int nc = wing3_section(w, 0, pxc, pyc);
    (void)wing3_section(w, 1, pxf, pyf);
    nring = nc;

    float halfflap = 0.5f * w->flap_frac * w->span;
    float zst[NST + 1];
    for (int t = 0; t <= NST; ++t) {
        float u = (float)t / NST;
        zst[t] = w->cz - 0.5f * w->span + u * w->span;
        int flapped = fabsf(zst[t] - w->cz) <= halfflap;
        const float *sx = flapped ? pxf : pxc;
        const float *sy = flapped ? pyf : pyc;
        for (int k = 0; k < nring; ++k) { rx[t][k] = sx[k]; ry[t][k] = sy[k]; }
    }

    float *o = out;
    int nv = 0;
    int cap = maxverts - 6; /* keep room for one quad = 2 tris = 6 verts */

    for (int t = 0; t < NST && nv < cap; ++t) {
        for (int k = 0; k < nring && nv < cap; ++k) {
            int k2 = (k + 1) % nring;
            float a[3] = { rx[t][k],   ry[t][k],   zst[t]   };
            float b[3] = { rx[t+1][k], ry[t+1][k], zst[t+1] };
            float c[3] = { rx[t+1][k2],ry[t+1][k2],zst[t+1] };
            float d[3] = { rx[t][k2],  ry[t][k2],  zst[t]   };
            emit_tri(&o, &nv, a, b, c);
            emit_tri(&o, &nv, a, c, d);
        }
    }
    /* tip caps: triangle fans about the section centroid */
    for (int side = 0; side < 2 && nv < cap; ++side) {
        int t = side == 0 ? 0 : NST;
        float czt[3] = { 0.f, 0.f, zst[t] };
        for (int k = 0; k < nring; ++k) { czt[0] += rx[t][k]; czt[1] += ry[t][k]; }
        czt[0] /= nring; czt[1] /= nring;
        for (int k = 0; k < nring && nv < cap; ++k) {
            int k2 = (k + 1) % nring;
            float a[3] = { rx[t][k],  ry[t][k],  zst[t] };
            float b[3] = { rx[t][k2], ry[t][k2], zst[t] };
            if (side == 0) emit_tri(&o, &nv, czt, b, a); /* outward -z */
            else           emit_tri(&o, &nv, czt, a, b); /* outward +z */
        }
    }
    return nv;
}

/* detailed visual mesh */
/* Everything below is render-only detail: the aerodynamic solid mask is
 * still the sealed plain-flap wing from wing3_rasterize. */

static float  *s_out, *s_mat;
static int     s_nv, s_max;
static float   s_cur[3];

static void em_setmat(float r, float g, float b)
{
    s_cur[0] = r; s_cur[1] = g; s_cur[2] = b;
}

static void em_tri(const float a[3], const float b[3], const float c[3])
{
    if (s_nv + 3 > s_max) return;
    float e1[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
    float e2[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
    float n[3]  = { e1[1]*e2[2]-e1[2]*e2[1],
                    e1[2]*e2[0]-e1[0]*e2[2],
                    e1[0]*e2[1]-e1[1]*e2[0] };
    float l = sqrtf(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (l > 1e-12f) { n[0] /= l; n[1] /= l; n[2] /= l; }
    const float *v[3] = { a, b, c };
    for (int k = 0; k < 3; ++k) {
        float *o = s_out + (size_t)s_nv * 6;
        o[0]=n[0]; o[1]=n[1]; o[2]=n[2];
        o[3]=v[k][0]; o[4]=v[k][1]; o[5]=v[k][2];
        float *m = s_mat + (size_t)s_nv * 3;
        m[0]=s_cur[0]; m[1]=s_cur[1]; m[2]=s_cur[2];
        ++s_nv;
    }
}

static void em_quad(const float a[3], const float b[3],
                    const float c[3], const float d[3])
{
    em_tri(a, b, c);
    em_tri(a, c, d);
}

/* section space (x/c, y/c) -> grid space, shared AoA/quarter-chord transform */
static void sec2grid(const Wing3 *w, float sx, float sy, float z, float out[3])
{
    float a = w->aoa_deg * PI_F / 180.f;
    float ca = cosf(a), sa = sinf(a);
    float X = sx * w->chord - 0.25f * w->chord;
    float Y = sy * w->chord;
    out[0] = w->cx + X * ca + Y * sa;
    out[1] = w->cy - X * sa + Y * ca;
    out[2] = z;
}

static void surf_pt(const Wing3 *w, float xc, int upper, float *sx, float *sy)
{
    float yc, dyc;
    camber_at(w, xc, &yc, &dyc);
    float yt = thickness_at(w, xc);
    float th = atanf(dyc), sn = sinf(th), cs = cosf(th);
    if (upper) { *sx = xc - yt * sn; *sy = yc + yt * cs; }
    else       { *sx = xc + yt * sn; *sy = yc - yt * cs; }
}

/* flap deflection in section space: rotation about the hinge plus a small
 * aft/down travel so a slot gap opens like a real slotted flap */
static void flap_pt(const Wing3 *w, float d, float sx, float sy,
                    float *ox, float *oy)
{
    float xh = w->flap_hinge, ych, dych;
    camber_at(w, xh, &ych, &dych);
    float k = fabsf(d) / (25.f * PI_F / 180.f);
    float cd = cosf(d), sd = sinf(d);
    float dx = sx - xh, dy = sy - ych;
    *ox = xh + dx * cd + dy * sd + 0.016f * k;
    *oy = ych - dx * sd + dy * cd - 0.010f * k;
}

/* extrude a constant section ring from z0 to z1 with optional end caps */
static void em_loft(const Wing3 *w, const float *rx, const float *ry, int n,
                    float z0, float z1, int cap0, int cap1)
{
    for (int k = 0; k < n; ++k) {
        int k2 = (k + 1) % n;
        float a[3], b[3], c[3], d[3];
        sec2grid(w, rx[k],  ry[k],  z0, a);
        sec2grid(w, rx[k],  ry[k],  z1, b);
        sec2grid(w, rx[k2], ry[k2], z1, c);
        sec2grid(w, rx[k2], ry[k2], z0, d);
        em_quad(a, b, c, d);
    }
    float cxm = 0.f, cym = 0.f;
    for (int k = 0; k < n; ++k) { cxm += rx[k]; cym += ry[k]; }
    cxm /= n; cym /= n;
    if (cap0) {
        float ce[3]; sec2grid(w, cxm, cym, z0, ce);
        for (int k = 0; k < n; ++k) {
            int k2 = (k + 1) % n;
            float a[3], b[3];
            sec2grid(w, rx[k], ry[k], z0, a);
            sec2grid(w, rx[k2], ry[k2], z0, b);
            em_tri(ce, b, a);
        }
    }
    if (cap1) {
        float ce[3]; sec2grid(w, cxm, cym, z1, ce);
        for (int k = 0; k < n; ++k) {
            int k2 = (k + 1) % n;
            float a[3], b[3];
            sec2grid(w, rx[k], ry[k], z1, a);
            sec2grid(w, rx[k2], ry[k2], z1, b);
            em_tri(ce, a, b);
        }
    }
}

/* tube between two grid-space points (actuator cylinders and rods) */
static void em_tube(const float p0[3], const float p1[3],
                    float r0, float r1, int nseg)
{
    float ax = p1[0]-p0[0], ay = p1[1]-p0[1], az = p1[2]-p0[2];
    float l = sqrtf(ax*ax + ay*ay + az*az);
    if (l < 1e-6f) return;
    ax /= l; ay /= l; az /= l;
    float ux, uy, uz;
    if (fabsf(ay) < 0.9f) { ux = -az; uy = 0.f; uz = ax; }
    else                  { ux = 0.f; uy = az; uz = -ay; }
    float ul = sqrtf(ux*ux + uy*uy + uz*uz);
    ux /= ul; uy /= ul; uz /= ul;
    float vx = ay*uz - az*uy, vy = az*ux - ax*uz, vz = ax*uy - ay*ux;
    for (int k = 0; k < nseg; ++k) {
        float t0 = 2.f*PI_F*k/nseg, t1 = 2.f*PI_F*(k+1)/nseg;
        float c0 = cosf(t0), s0 = sinf(t0), c1 = cosf(t1), s1 = sinf(t1);
        float a[3] = { p0[0]+r0*(c0*ux+s0*vx), p0[1]+r0*(c0*uy+s0*vy), p0[2]+r0*(c0*uz+s0*vz) };
        float b[3] = { p1[0]+r1*(c0*ux+s0*vx), p1[1]+r1*(c0*uy+s0*vy), p1[2]+r1*(c0*uz+s0*vz) };
        float c[3] = { p1[0]+r1*(c1*ux+s1*vx), p1[1]+r1*(c1*uy+s1*vy), p1[2]+r1*(c1*uz+s1*vz) };
        float d[3] = { p0[0]+r0*(c1*ux+s1*vx), p0[1]+r0*(c1*uy+s1*vy), p0[2]+r0*(c1*uz+s1*vz) };
        em_quad(a, b, c, d);
        em_tri(p1, c, b); /* end disks */
        em_tri(p0, a, d);
    }
}

/* elongated fairing pod under the flap tracks. The part of the canoe aft
 * of the hinge articulates down with the flap, like the real mechanism --
 * otherwise the tails would stick straight out into the air. */
static void pod_pt(const Wing3 *w, float d, float x, float y, float z,
                   float dz, float out[3])
{
    if (x > w->flap_hinge && d != 0.f) {
        float fx, fy;
        flap_pt(w, d, x, y, &fx, &fy);
        sec2grid(w, fx, fy, z, out);
    } else {
        sec2grid(w, x, y, z, out);
    }
    out[2] += dz;
}

static void em_pod(const Wing3 *w, float d, float xc0, float xc1, float ybase,
                   float z, float ry, float rz)
{
    enum { NU = 12, NR = 10 };
    for (int i = 0; i < NU; ++i) {
        float u0 = (float)i / NU, u1 = (float)(i + 1) / NU;
        float x0 = xc0 + (xc1 - xc0) * u0, x1 = xc0 + (xc1 - xc0) * u1;
        float w0 = sinf(PI_F * u0), w1 = sinf(PI_F * u1); /* pointy ends */
        for (int k = 0; k < NR; ++k) {
            float t0 = PI_F * k / NR, t1 = PI_F * (k + 1) / NR;
            float a[3], b[3], c[3], d4[3];
            pod_pt(w, d, x0, ybase - ry*w0*sinf(t0), z, rz*w0*cosf(t0), a);
            pod_pt(w, d, x1, ybase - ry*w1*sinf(t0), z, rz*w1*cosf(t0), b);
            pod_pt(w, d, x1, ybase - ry*w1*sinf(t1), z, rz*w1*cosf(t1), c);
            pod_pt(w, d, x0, ybase - ry*w0*sinf(t1), z, rz*w0*cosf(t1), d4);
            em_quad(a, b, c, d4);
        }
    }
}

int wing3_mesh_detailed(const Wing3 *w, float *out, float *matcol, int maxverts)
{
    s_out = out; s_mat = matcol; s_nv = 0; s_max = maxverts;

    const float d   = w->flap_deg * PI_F / 180.f;
    const float xh  = w->flap_hinge;
    const float zf0 = w->cz - 0.5f * w->flap_frac * w->span;
    const float zf1 = w->cz + 0.5f * w->flap_frac * w->span;
    const float zt0 = w->cz - 0.5f * w->span;
    const float zt1 = w->cz + 0.5f * w->span;
    enum { NS2 = 46 };
    float rx[2 * NS2 + 8], ry[2 * NS2 + 8];

    /* full section ring: outboard panels */
    int nfull = 0;
    for (int k = NS2 - 1; k >= 0; --k) {
        float xc = 0.5f * (1.f - cosf(PI_F * k / (NS2 - 1)));
        surf_pt(w, xc, 1, &rx[nfull], &ry[nfull]); ++nfull;
    }
    for (int k = 1; k < NS2; ++k) {
        float xc = 0.5f * (1.f - cosf(PI_F * k / (NS2 - 1)));
        surf_pt(w, xc, 0, &rx[nfull], &ry[nfull]); ++nfull;
    }
    em_setmat(0.84f, 0.85f, 0.88f);
    em_loft(w, rx, ry, nfull, zt0, zf0 - 0.15f, 1, 1);
    em_loft(w, rx, ry, nfull, zf1 + 0.15f, zt1, 1, 1);

    /* main element with cove cut (flap bay) */
    const float lip = xh + 0.06f; /* upper skin overhangs the flap  */
    const float low = xh - 0.03f; /* lower skin stops short         */
    int nmain = 0;
    for (int k = NS2 - 1; k >= 0; --k) { /* upper: lip -> LE */
        float xc = lip * 0.5f * (1.f - cosf(PI_F * k / (NS2 - 1)));
        surf_pt(w, xc, 1, &rx[nmain], &ry[nmain]); ++nmain;
    }
    for (int k = 1; k < NS2; ++k) {      /* lower: LE -> low */
        float xc = low * 0.5f * (1.f - cosf(PI_F * k / (NS2 - 1)));
        surf_pt(w, xc, 0, &rx[nmain], &ry[nmain]); ++nmain;
    }
    { /* cove: concave return from the lower lip up into the bay */
        float ux2, uy2, lx2, ly2, ych, dych;
        surf_pt(w, lip, 1, &ux2, &uy2);
        surf_pt(w, low, 0, &lx2, &ly2);
        camber_at(w, xh, &ych, &dych);
        rx[nmain] = low + 0.01f;  ry[nmain] = ly2 + 0.35f * (ych - ly2); ++nmain;
        rx[nmain] = xh - 0.005f;  ry[nmain] = ych + 0.55f * (uy2 - ych); ++nmain;
        rx[nmain] = lip - 0.004f; ry[nmain] = uy2 - 0.012f;              ++nmain;
    }
    em_setmat(0.84f, 0.85f, 0.88f);
    em_loft(w, rx, ry, nmain, zf0, zf1, 1, 1);

    /* flap body (deflected, slot gap opens with deflection) */
    int nflap = 0;
    const float fu = xh + 0.018f, fl = xh - 0.012f;
    for (int k = NS2 - 1; k >= 0; --k) { /* upper: TE -> flap LE */
        float xc = fu + (1.f - fu) * 0.5f * (1.f - cosf(PI_F * k / (NS2 - 1)));
        float sx, sy; surf_pt(w, xc, 1, &sx, &sy);
        flap_pt(w, d, sx, sy, &rx[nflap], &ry[nflap]); ++nflap;
    }
    { /* rounded flap leading edge */
        float u0x, u0y, l0x, l0y;
        surf_pt(w, fu, 1, &u0x, &u0y);
        surf_pt(w, fl, 0, &l0x, &l0y);
        for (int k = 1; k <= 4; ++k) {
            float t = (float)k / 5.f;
            float mx = u0x + (l0x - u0x) * t - 0.016f * sinf(PI_F * t);
            float my = u0y + (l0y - u0y) * t;
            flap_pt(w, d, mx, my, &rx[nflap], &ry[nflap]); ++nflap;
        }
    }
    for (int k = 1; k < NS2; ++k) {      /* lower: flap LE -> TE */
        float xc = fl + (1.f - fl) * 0.5f * (1.f - cosf(PI_F * k / (NS2 - 1)));
        float sx, sy; surf_pt(w, xc, 0, &sx, &sy);
        flap_pt(w, d, sx, sy, &rx[nflap], &ry[nflap]); ++nflap;
    }
    em_setmat(0.78f, 0.80f, 0.84f);
    em_loft(w, rx, ry, nflap, zf0 + 0.3f, zf1 - 0.3f, 1, 1);

    int nv_paint = s_nv; /* panel lines apply to the skin only */

    /* flap track fairings + hydraulic actuators */
    const float zpods[3] = { w->cz - 0.30f * w->flap_frac * w->span,
                             w->cz,
                             w->cz + 0.30f * w->flap_frac * w->span };
    for (int p = 0; p < 3; ++p) {
        float ylow, sxx;
        surf_pt(w, 0.62f, 0, &sxx, &ylow);
        em_setmat(0.66f, 0.67f, 0.72f);
        em_pod(w, d, 0.45f, 1.04f, ylow + 0.004f, zpods[p], 0.030f, 2.0f);

        /* actuator: barrel fixed under the main element, polished rod out
         * to the flap -- it visibly extends as the flap goes down */
        float ax, ay;
        surf_pt(w, 0.52f, 0, &ax, &ay);
        float A[3]; sec2grid(w, ax, ay - 0.012f, zpods[p], A);
        float bx, by, sx2, sy2;
        surf_pt(w, xh + 0.16f, 0, &sx2, &sy2);
        flap_pt(w, d, sx2, sy2 - 0.006f, &bx, &by);
        float B[3]; sec2grid(w, bx, by, zpods[p], B);
        float M[3] = { A[0] + (B[0]-A[0]) * 0.55f,
                       A[1] + (B[1]-A[1]) * 0.55f,
                       A[2] + (B[2]-A[2]) * 0.55f };
        em_setmat(0.30f, 0.31f, 0.34f);
        em_tube(A, M, 0.85f, 0.85f, 10);
        em_setmat(0.92f, 0.94f, 0.97f);
        em_tube(M, B, 0.38f, 0.38f, 8);
    }

    /* vortex generators: alternating angled fins at 20% chord */
    em_setmat(0.55f, 0.56f, 0.60f);
    {
        float sx, sy; surf_pt(w, 0.20f, 1, &sx, &sy);
        int nvg = 12;
        for (int k = 0; k < nvg; ++k) {
            float z = zf0 + (zf1 - zf0) * (k + 0.5f) / nvg;
            float yaw = (k & 1) ? 0.30f : -0.30f;
            float a[3], b[3], c[3];
            sec2grid(w, sx, sy, z, a);
            sec2grid(w, sx + 0.055f, sy + 0.003f, z, b);
            b[2] += 2.2f * sinf(yaw);
            sec2grid(w, sx + 0.050f, sy + 0.030f, z, c);
            em_tri(a, b, c);
            em_tri(a, c, b); /* visible from both sides */
        }
    }

    /* baked panel lines on the skin */
    for (int v = 0; v < nv_paint; ++v) {
        float *pos = s_out + (size_t)v * 6 + 3;
        float *m   = s_mat + (size_t)v * 3;
        if (fmodf(pos[2], 8.f) < 0.55f || fmodf(pos[0], 11.f) < 0.55f) {
            m[0] *= 0.86f; m[1] *= 0.86f; m[2] *= 0.86f;
        }
    }
    return s_nv;
}

void wing3_lights(const Wing3 *w, float pos[4][3])
{
    float sx, sy;
    surf_pt(w, 0.10f, 1, &sx, &sy);
    sec2grid(w, sx, sy, w->cz - 0.5f * w->span - 0.6f, pos[0]); /* port red   */
    sec2grid(w, sx, sy, w->cz + 0.5f * w->span + 0.6f, pos[1]); /* stbd green */
    surf_pt(w, 0.92f, 1, &sx, &sy);
    sec2grid(w, sx, sy + 0.01f, w->cz - 0.5f * w->span - 0.4f, pos[2]);
    sec2grid(w, sx, sy + 0.01f, w->cz + 0.5f * w->span + 0.4f, pos[3]);
}
