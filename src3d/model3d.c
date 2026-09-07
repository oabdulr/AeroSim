#include "model3d.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_TRI 600000

static float *g_raw;      /* as loaded: 9 floats per triangle              */
static float *g_fit;      /* axis-permuted, centered, scaled to grid units */
static int    g_ntri;
static int    g_axmode;   /* 0: xyz, 1: x z -y (z-up files), 2: z y -x     */
static char   g_name[128];
static double g_refS = 1.0;
static float  g_cx, g_cy, g_cz; /* placement in the grid */

bool model3_active(void) { return g_ntri > 0; }
const char *model3_name(void) { return g_name; }
double model3_refS(void) { return g_refS; }
int model3_tris(void) { return g_ntri; }
void model3_clear(void) { g_ntri = 0; }

/* parsing */

static bool parse_stl_binary(const uint8_t *buf, long len)
{
    if (len < 84) return false;
    uint32_t n;
    memcpy(&n, buf + 80, 4);
    if (n == 0 || (long)(84 + (long long)n * 50) != len) return false;
    if (n > MAX_TRI) n = MAX_TRI;
    g_ntri = (int)n;
    for (uint32_t t = 0; t < n; ++t) {
        const uint8_t *rec = buf + 84 + (size_t)t * 50;
        memcpy(&g_raw[t * 9], rec + 12, 36); /* skip the stored normal */
    }
    return true;
}

static bool parse_stl_ascii(const char *txt)
{
    const char *p = txt;
    int nv = 0;
    g_ntri = 0;
    while ((p = strstr(p, "vertex")) != NULL) {
        float x, y, z;
        if (sscanf(p + 6, "%f %f %f", &x, &y, &z) == 3) {
            if (g_ntri >= MAX_TRI) break;
            float *dst = &g_raw[g_ntri * 9 + (nv % 3) * 3];
            dst[0] = x; dst[1] = y; dst[2] = z;
            ++nv;
            if (nv % 3 == 0) ++g_ntri;
        }
        p += 6;
    }
    return g_ntri > 0;
}

static bool parse_obj(const char *txt)
{
    /* Pass 1: count vertices. */
    int nverts = 0;
    for (const char *p = txt; p; p = strchr(p, '\n')) {
        while (*p == '\n' || *p == '\r') ++p;
        if (p[0] == 'v' && (p[1] == ' ' || p[1] == '\t')) ++nverts;
        if (!*p) break;
    }
    if (nverts == 0) return false;
    float *vx = (float *)malloc(sizeof(float) * 3 * (size_t)nverts);
    if (!vx) return false;

    int vi = 0;
    g_ntri = 0;
    for (const char *p = txt; p; ) {
        while (*p == '\n' || *p == '\r') ++p;
        if (!*p) break;
        if (p[0] == 'v' && (p[1] == ' ' || p[1] == '\t')) {
            if (vi < nverts &&
                sscanf(p + 1, "%f %f %f", &vx[vi*3], &vx[vi*3+1], &vx[vi*3+2]) == 3)
                ++vi;
        } else if (p[0] == 'f' && (p[1] == ' ' || p[1] == '\t')) {
            /* triangulate the polygon as a fan */
            int idx[64], m = 0;
            const char *q = p + 1;
            while (*q && *q != '\n' && m < 64) {
                while (*q == ' ' || *q == '\t') ++q;
                if (*q == '\n' || *q == '\r' || !*q) break;
                int v = atoi(q);
                if (v < 0) v = vi + v + 1;
                if (v >= 1 && v <= vi) idx[m++] = v - 1;
                while (*q && *q != ' ' && *q != '\t' && *q != '\n') ++q;
            }
            for (int k = 2; k < m && g_ntri < MAX_TRI; ++k) {
                float *dst = &g_raw[g_ntri * 9];
                memcpy(dst,     &vx[idx[0]     * 3], 12);
                memcpy(dst + 3, &vx[idx[k - 1] * 3], 12);
                memcpy(dst + 6, &vx[idx[k]     * 3], 12);
                ++g_ntri;
            }
        }
        p = strchr(p, '\n');
    }
    free(vx);
    return g_ntri > 0;
}

/* Axis permutation, then center + uniform scale to fit the tunnel. */
static void refit(void)
{
    for (int t = 0; t < g_ntri; ++t) {
        for (int v = 0; v < 3; ++v) {
            const float *s = &g_raw[t * 9 + v * 3];
            float *d = &g_fit[t * 9 + v * 3];
            switch (g_axmode) {
            case 0: d[0] = s[0]; d[1] = s[1]; d[2] = s[2]; break;
            case 1: d[0] = s[0]; d[1] = s[2]; d[2] = -s[1]; break; /* z-up */
            default: d[0] = s[2]; d[1] = s[1]; d[2] = -s[0]; break;
            }
        }
    }
    float lo[3] = {1e30f, 1e30f, 1e30f}, hi[3] = {-1e30f, -1e30f, -1e30f};
    for (int k = 0; k < g_ntri * 3; ++k)
        for (int c = 0; c < 3; ++c) {
            float v = g_fit[k * 3 + c];
            if (v < lo[c]) lo[c] = v;
            if (v > hi[c]) hi[c] = v;
        }
    float ex = hi[0] - lo[0], ey = hi[1] - lo[1], ez = hi[2] - lo[2];
    if (ex < 1e-9f) ex = 1e-9f;
    if (ey < 1e-9f) ey = 1e-9f;
    if (ez < 1e-9f) ez = 1e-9f;
    /* keep the model inside the live region, clear of the absorbing layers */
    float sc = 90.f / ex;
    if (44.f / ey < sc) sc = 44.f / ey;
    if (60.f / ez < sc) sc = 60.f / ez;
    /* blockage cap: frontal area <= ~6% of the tunnel cross-section.
     * Bigger than that and the squeezed flow inflates the measured CD
     * (classic wind-tunnel blockage error) -- a face-on cube read 2.5
     * instead of ~1.1 before this cap. */
    {
        float fmax = 0.06f * (float)NY3 * (float)NZ3;
        float frontal = (ey * sc) * (ez * sc);
        if (frontal > fmax) sc *= sqrtf(fmax / frontal);
    }
    float cx = 0.5f * (lo[0] + hi[0]);
    float cy = 0.5f * (lo[1] + hi[1]);
    float cz = 0.5f * (lo[2] + hi[2]);
    for (int k = 0; k < g_ntri * 3; ++k) {
        g_fit[k * 3 + 0] = (g_fit[k * 3 + 0] - cx) * sc;
        g_fit[k * 3 + 1] = (g_fit[k * 3 + 1] - cy) * sc;
        g_fit[k * 3 + 2] = (g_fit[k * 3 + 2] - cz) * sc;
    }
    g_refS = (double)(ey * sc) * (double)(ez * sc); /* frontal-area proxy */
    g_cx = NX3 * 0.34f;
    g_cy = NY3 * 0.5f;
    g_cz = NZ3 * 0.5f;
}

bool model3_load(const char *path, char *err, int errlen)
{
    if (!g_raw) g_raw = (float *)malloc(sizeof(float) * 9 * MAX_TRI);
    if (!g_fit) g_fit = (float *)malloc(sizeof(float) * 9 * MAX_TRI);
    if (!g_raw || !g_fit) { snprintf(err, errlen, "out of memory"); return false; }

    FILE *f = fopen(path, "rb");
    if (!f) { snprintf(err, errlen, "cannot open file"); return false; }
    fseek(f, 0, SEEK_END);
    long len = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (len <= 0 || len > 400 * 1024 * 1024) {
        fclose(f);
        snprintf(err, errlen, "file empty or over 400 MB");
        return false;
    }
    uint8_t *buf = (uint8_t *)malloc((size_t)len + 1);
    if (!buf) { fclose(f); snprintf(err, errlen, "out of memory"); return false; }
    size_t rd = fread(buf, 1, (size_t)len, f);
    fclose(f);
    buf[rd] = 0;

    g_ntri = 0;
    const char *dot = strrchr(path, '.');
    bool ok = false;
    if (dot && (_stricmp(dot, ".obj") == 0)) {
        ok = parse_obj((const char *)buf);
    } else {
        ok = parse_stl_binary(buf, (long)rd);
        if (!ok) ok = parse_stl_ascii((const char *)buf);
    }
    free(buf);

    if (!ok || g_ntri == 0) {
        g_ntri = 0;
        snprintf(err, errlen, "could not parse (supported: binary/ASCII STL, OBJ)");
        return false;
    }

    const char *base = strrchr(path, '\\');
    if (!base) base = strrchr(path, '/');
    snprintf(g_name, sizeof g_name, "%s", base ? base + 1 : path);
    g_axmode = 0;
    refit();
    return true;
}

void model3_cycle_axes(void)
{
    if (!model3_active()) return;
    g_axmode = (g_axmode + 1) % 3;
    refit();
}

/* placement transform */
/* yaw about the vertical axis first, then nose-up pitch (clockwise in the
 * x-y plane, y-up), then translate into the grid */
static void xform(const float *v, float ca, float sa, float cy2, float sy2,
                  float *out)
{
    float x = v[0] * cy2 + v[2] * sy2;
    float z = -v[0] * sy2 + v[2] * cy2;
    out[0] = g_cx + x * ca + v[1] * sa;
    out[1] = g_cy - x * sa + v[1] * ca;
    out[2] = g_cz + z;
}

/* voxelization: x-direction ray parity fill */
#define MAXCROSS 32

void model3_rasterize(Sim3 *s, float pitch_deg, float yaw_deg)
{
    static uint8_t *news = NULL;
    static float   *xs = NULL;   /* [NY3*NZ3][MAXCROSS] crossing coordinates */
    static uint8_t *xn = NULL;
    if (!news) news = (uint8_t *)malloc(NCELLS3);
    if (!xs)   xs = (float *)malloc(sizeof(float) * MAXCROSS * NY3 * NZ3);
    if (!xn)   xn = (uint8_t *)malloc((size_t)NY3 * NZ3);
    memset(news, 0, NCELLS3);
    memset(xn, 0, (size_t)NY3 * NZ3);

    float a = pitch_deg * 3.14159265f / 180.f;
    float ca = cosf(a), sa = sinf(a);
    float ya = yaw_deg * 3.14159265f / 180.f;
    float cy2 = cosf(ya), sy2 = sinf(ya);

    /* tiny irrational ray offsets dodge exact edge/vertex hits, whose double
     * counting would flip the fill parity */
    const float oy = 0.00731f, oz = 0.00577f;

    for (int t = 0; t < g_ntri; ++t) {
        float A[3], B[3], C[3];
        xform(&g_fit[t * 9 + 0], ca, sa, cy2, sy2, A);
        xform(&g_fit[t * 9 + 3], ca, sa, cy2, sy2, B);
        xform(&g_fit[t * 9 + 6], ca, sa, cy2, sy2, C);

        float ylo = A[1], yhi = A[1], zlo = A[2], zhi = A[2];
        if (B[1] < ylo) ylo = B[1]; if (B[1] > yhi) yhi = B[1];
        if (C[1] < ylo) ylo = C[1]; if (C[1] > yhi) yhi = C[1];
        if (B[2] < zlo) zlo = B[2]; if (B[2] > zhi) zhi = B[2];
        if (C[2] < zlo) zlo = C[2]; if (C[2] > zhi) zhi = C[2];

        int y0 = (int)ceilf(ylo - oy), y1 = (int)floorf(yhi - oy);
        int z0 = (int)ceilf(zlo - oz), z1 = (int)floorf(zhi - oz);
        if (y0 < 1) y0 = 1; if (y1 > NY3 - 2) y1 = NY3 - 2;
        if (z0 < 1) z0 = 1; if (z1 > NZ3 - 2) z1 = NZ3 - 2;

        float d = (B[1] - A[1]) * (C[2] - A[2]) - (B[2] - A[2]) * (C[1] - A[1]);
        if (fabsf(d) < 1e-12f) continue; /* edge-on to the ray direction */
        float inv = 1.f / d;

        for (int z = z0; z <= z1; ++z) {
            float pz = (float)z + oz;
            for (int y = y0; y <= y1; ++y) {
                float py = (float)y + oy;
                float u = ((py - A[1]) * (C[2] - A[2]) - (pz - A[2]) * (C[1] - A[1])) * inv;
                float v = ((B[1] - A[1]) * (pz - A[2]) - (B[2] - A[2]) * (py - A[1])) * inv;
                if (u < 0.f || v < 0.f || u + v > 1.f) continue;
                float x = A[0] + u * (B[0] - A[0]) + v * (C[0] - A[0]);
                int r = z * NY3 + y;
                if (xn[r] < MAXCROSS) xs[r * MAXCROSS + xn[r]++] = x;
            }
        }
    }

    for (int z = 1; z < NZ3 - 1; ++z) {
        for (int y = 1; y < NY3 - 1; ++y) {
            int r = z * NY3 + y;
            int m = xn[r];
            if (m < 2) continue;
            float *c = &xs[r * MAXCROSS];
            for (int i2 = 1; i2 < m; ++i2) { /* insertion sort */
                float v = c[i2]; int b = i2 - 1;
                while (b >= 0 && c[b] > v) { c[b + 1] = c[b]; --b; }
                c[b + 1] = v;
            }
            for (int k = 0; k + 1 < m; k += 2) {
                int xa = (int)ceilf(c[k]);
                int xb = (int)floorf(c[k + 1]);
                if (xa < 1) xa = 1;
                if (xb > NX3 - 2) xb = NX3 - 2;
                for (int x = xa; x <= xb; ++x) news[IDX3(x, y, z)] = 1;
            }
        }
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

int model3_mesh(float *out, int maxverts, float pitch_deg, float yaw_deg)
{
    float a = pitch_deg * 3.14159265f / 180.f;
    float ca = cosf(a), sa = sinf(a);
    float ya = yaw_deg * 3.14159265f / 180.f;
    float cy2 = cosf(ya), sy2 = sinf(ya);
    int nv = 0;
    float *o = out;
    for (int t = 0; t < g_ntri && nv + 3 <= maxverts; ++t) {
        float A[3], B[3], C[3];
        xform(&g_fit[t * 9 + 0], ca, sa, cy2, sy2, A);
        xform(&g_fit[t * 9 + 3], ca, sa, cy2, sy2, B);
        xform(&g_fit[t * 9 + 6], ca, sa, cy2, sy2, C);
        float e1[3] = {B[0]-A[0], B[1]-A[1], B[2]-A[2]};
        float e2[3] = {C[0]-A[0], C[1]-A[1], C[2]-A[2]};
        float n[3] = {e1[1]*e2[2]-e1[2]*e2[1],
                      e1[2]*e2[0]-e1[0]*e2[2],
                      e1[0]*e2[1]-e1[1]*e2[0]};
        float l = sqrtf(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if (l > 1e-12f) { n[0] /= l; n[1] /= l; n[2] /= l; }
        const float *V[3] = {A, B, C};
        for (int k = 0; k < 3; ++k) {
            o[0]=n[0]; o[1]=n[1]; o[2]=n[2];
            o[3]=V[k][0]; o[4]=V[k][1]; o[5]=V[k][2];
            o += 6; ++nv;
        }
    }
    return nv;
}
