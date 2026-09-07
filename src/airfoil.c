#include "airfoil.h"
#include <math.h>
#include <string.h>

#define PI_F 3.14159265358979f

/* NACA 4-digit mean camber line and its slope at chordwise position xc. */
static void camber_at(const Airfoil *af, float xc, float *yc, float *dyc)
{
    float m = af->camber, p = af->camber_pos;
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

/* NACA 4-digit half thickness (closed-trailing-edge coefficient -0.1036). */
static float thickness_at(const Airfoil *af, float xc)
{
    float t = af->thickness;
    return 5.f * t * (0.2969f * sqrtf(xc) - 0.1260f * xc - 0.3516f * xc * xc
                      + 0.2843f * xc * xc * xc - 0.1036f * xc * xc * xc * xc);
}

void airfoil_build(Airfoil *af)
{
    enum { NS = 160 }; /* samples per surface; 2*NS-1 <= AF_MAX_PTS */
    float sxu[NS], syu[NS], sxl[NS], syl[NS];

    /* Plain flap: everything aft of the hinge rotates about the hinge point
     * on the camber line. Positive deflection = trailing edge down. */
    float d  = af->flap_deg * PI_F / 180.f;
    float xh = af->flap_hinge;
    float ych, dych;
    camber_at(af, xh, &ych, &dych);
    float cd = cosf(d), sd = sinf(d);

    for (int k = 0; k < NS; ++k) {
        float sf = (float)k / (float)(NS - 1);
        float xc = 0.5f * (1.f - cosf(PI_F * sf)); /* cosine clustering */
        float yc, dyc;
        camber_at(af, xc, &yc, &dyc);
        float yt = thickness_at(af, xc);
        float th = atanf(dyc);
        float sn = sinf(th), cs = cosf(th);

        /* Thickness applied perpendicular to the camber line. */
        float xu = xc - yt * sn, yu = yc + yt * cs;
        float xl = xc + yt * sn, yl = yc - yt * cs;

        if (xu > xh) {
            float dx = xu - xh, dy = yu - ych;
            xu = xh + dx * cd + dy * sd;
            yu = ych - dx * sd + dy * cd;
        }
        if (xl > xh) {
            float dx = xl - xh, dy = yl - ych;
            xl = xh + dx * cd + dy * sd;
            yl = ych - dx * sd + dy * cd;
        }
        sxu[k] = xu; syu[k] = yu;
        sxl[k] = xl; syl[k] = yl;
    }

    /* Closed outline: upper surface TE -> LE, then lower surface LE -> TE. */
    int n = 0;
    for (int k = NS - 1; k >= 0; --k) { af->px[n] = sxu[k]; af->py[n] = syu[k]; ++n; }
    for (int k = 1; k < NS; ++k)      { af->px[n] = sxl[k]; af->py[n] = syl[k]; ++n; }
    af->npts = n;

    /* Scale to chord, rotate by AoA about the quarter-chord point, translate.
     * Nose-up is a clockwise rotation in y-up grid coordinates. */
    float a  = af->aoa_deg * PI_F / 180.f;
    float ca = cosf(a), sa = sinf(a);
    for (int k = 0; k < n; ++k) {
        float X = af->px[k] * af->chord - 0.25f * af->chord;
        float Y = af->py[k] * af->chord;
        af->px[k] = af->cx + X * ca + Y * sa;
        af->py[k] = af->cy - X * sa + Y * ca;
    }
}

/* Scanline fill with the nonzero winding rule. Winding (rather than even-odd)
 * matters here: a deflected plain flap can make the outline self-overlap
 * slightly near the hinge, and even-odd would punch a leak hole there. */
void airfoil_rasterize(const Airfoil *af, Sim *s)
{
    static uint8_t news[NCELLS];
    memset(news, 0, NCELLS);

    const float *px = af->px, *py = af->py;
    const int    n  = af->npts;

    float ymin = 1e9f, ymax = -1e9f;
    for (int k = 0; k < n; ++k) {
        if (py[k] < ymin) ymin = py[k];
        if (py[k] > ymax) ymax = py[k];
    }
    int y0 = (int)floorf(ymin) - 1; if (y0 < 1) y0 = 1;
    int y1 = (int)ceilf(ymax) + 1;  if (y1 > NY - 2) y1 = NY - 2;

    for (int y = y0; y <= y1; ++y) {
        float xs[64];
        int   dir[64];
        int   m  = 0;
        float fy = (float)y;

        for (int k = 0; k < n; ++k) {
            float ax = px[k], ay = py[k];
            float bx = px[(k + 1) % n], by = py[(k + 1) % n];
            if ((ay > fy) != (by > fy)) {
                float t = (fy - ay) / (by - ay);
                if (m < 64) {
                    xs[m]  = ax + t * (bx - ax);
                    dir[m] = (by > ay) ? 1 : -1;
                    ++m;
                }
            }
        }
        /* Insertion sort of crossings by x (m is tiny). */
        for (int a = 1; a < m; ++a) {
            float vx = xs[a]; int vd = dir[a];
            int b = a - 1;
            while (b >= 0 && xs[b] > vx) {
                xs[b + 1] = xs[b]; dir[b + 1] = dir[b]; --b;
            }
            xs[b + 1] = vx; dir[b + 1] = vd;
        }
        /* Fill spans where the winding number is nonzero. */
        int wind = 0;
        for (int k = 0; k + 1 < m; ++k) {
            wind += dir[k];
            if (wind != 0) {
                int xa = (int)ceilf(xs[k]);
                int xb = (int)floorf(xs[k + 1]);
                if (xa < 1) xa = 1;
                if (xb > NX - 2) xb = NX - 2;
                for (int x = xa; x <= xb; ++x) news[y * NX + x] = 1;
            }
        }
    }

    /* Cells that just switched from solid to fluid need valid populations;
     * cells that became solid get neutral macro fields so the vorticity
     * stencil and arrows never see stale velocities next to the surface. */
    for (int i = 0; i < NCELLS; ++i) {
        if (s->solid[i] && !news[i]) {
            sim_refresh_cell(s, i);
        } else if (!s->solid[i] && news[i]) {
            s->rho[i] = 1.f;
            s->ux[i]  = 0.f;
            s->uy[i]  = 0.f;
        }
    }

    memcpy(s->solid, news, NCELLS);
    sim_update_careful_mask(s);
}
