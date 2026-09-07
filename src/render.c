#include "render.h"
#include <math.h>
#include <stdlib.h>

static inline float clamp01(float t) { return t < 0.f ? 0.f : (t > 1.f ? 1.f : t); }

static inline uint32_t pack(float r, float g, float b)
{
    int R = (int)(clamp01(r) * 255.f + 0.5f);
    int G = (int)(clamp01(g) * 255.f + 0.5f);
    int B = (int)(clamp01(b) * 255.f + 0.5f);
    return ((uint32_t)R << 16) | ((uint32_t)G << 8) | (uint32_t)B;
}

/* Rainbow map (CFD "jet"): dark blue -> cyan -> green -> yellow -> red. */
static uint32_t cmap_jet(float t)
{
    t = clamp01(t);
    float r = clamp01(1.5f - fabsf(4.f * t - 3.f));
    float g = clamp01(1.5f - fabsf(4.f * t - 2.f));
    float b = clamp01(1.5f - fabsf(4.f * t - 1.f));
    /* darken the extremes slightly so 0 and 1 stay readable */
    if (t < 0.125f) b = 0.55f + 3.6f * t;
    if (t > 0.875f) r = 1.f - (t - 0.875f) * 1.6f;
    return pack(r, g, b);
}

/* Diverging blue -> white -> red (for signed fields like vorticity). */
static uint32_t cmap_div(float t)
{
    t = clamp01(t);
    float r, g, b;
    if (t < 0.5f) {
        float u = t * 2.f;
        r = 0.12f + (0.97f - 0.12f) * u;
        g = 0.22f + (0.97f - 0.22f) * u;
        b = 0.83f + (0.97f - 0.83f) * u;
    } else {
        float u = (t - 0.5f) * 2.f;
        r = 0.97f + (0.79f - 0.97f) * u;
        g = 0.97f + (0.09f - 0.97f) * u;
        b = 0.97f + (0.12f - 0.97f) * u;
    }
    return pack(r, g, b);
}

const char *vis_name(int mode)
{
    switch (mode) {
    case VIS_PRESSURE:  return "Pressure (Cp)";
    case VIS_SPEED:     return "Velocity magnitude";
    case VIS_VORTICITY: return "Vorticity";
    }
    return "?";
}

const char *vis_range(int mode)
{
    switch (mode) {
    case VIS_PRESSURE:  return "Cp  -2.0 (low)  to  +1.0 (high / stagnation)";
    case VIS_SPEED:     return "|u| / U0   0.0  to  1.6";
    case VIS_VORTICITY: return "clockwise (blue)  to  counter-clockwise (red)";
    }
    return "";
}

void render_field(const Sim *s, uint32_t *pix, int mode)
{
    const float u0 = s->u0;
    int i, y; /* MSVC's OpenMP wants parallel loop indices declared here */

    if (mode == VIS_PRESSURE) {
        /* Cp = (p - p_inf) / (0.5 rho0 u0^2), with p = cs^2 rho. Only
         * pressure differences are physical in LBM, so reference to the
         * undisturbed freestream: an upstream strip inside the absorbing
         * layer, where the flow is anchored to far-field conditions. */
        double msum = 0.0;
        int    mcount = 0;
        for (int yy = NY / 2 - 40; yy < NY / 2 + 40; ++yy) {
            for (int xx = 4; xx < 12; ++xx) {
                msum += s->rho[yy * NX + xx];
                ++mcount;
            }
        }
        const float rho_ref = (float)(msum / (double)mcount);
        const float inv_qd  = 1.f / (0.5f * u0 * u0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (i = 0; i < NCELLS; ++i) {
            float cp = ((s->rho[i] - rho_ref) / 3.f) * inv_qd;
            pix[i] = cmap_jet((cp + 2.f) / 3.f); /* map [-2,+1] -> [0,1] */
        }
    } else if (mode == VIS_SPEED) {
        const float inv = 1.f / (1.6f * u0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (i = 0; i < NCELLS; ++i) {
            float vx = s->ux[i], vy = s->uy[i];
            pix[i] = cmap_jet(sqrtf(vx * vx + vy * vy) * inv);
        }
    } else {
        /* Vorticity by central differences, scaled to ~ +-25 U0/chord. */
        const float wscale = 25.f * u0 / s->chord;
        const float inv2w  = 0.5f / wscale;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (y = 0; y < NY; ++y) {
            for (int x = 0; x < NX; ++x) {
                int i = y * NX + x;
                if (x == 0 || x == NX - 1 || y == 0 || y == NY - 1 || s->solid[i]) {
                    pix[i] = cmap_div(0.5f);
                    continue;
                }
                float w = 0.5f * (s->uy[i + 1] - s->uy[i - 1])
                        - 0.5f * (s->ux[i + NX] - s->ux[i - NX]);
                pix[i] = cmap_div(0.5f + w * inv2w);
            }
        }
    }
}

static void draw_line_f(uint32_t *pix, float x0, float y0, float x1, float y1, uint32_t c)
{
    float dx = x1 - x0, dy = y1 - y0;
    int   ns = (int)fmaxf(fabsf(dx), fabsf(dy)) + 1;
    for (int k = 0; k <= ns; ++k) {
        float t = (float)k / (float)ns;
        int   x = (int)(x0 + dx * t + 0.5f);
        int   y = (int)(y0 + dy * t + 0.5f);
        if (x >= 0 && x < NX && y >= 0 && y < NY) pix[y * NX + x] = c;
    }
}

void render_airfoil(const Sim *s, const Airfoil *af, uint32_t *pix)
{
    const uint32_t body = 0x00E9E9EE, edge = 0x0020222A;
    for (int i = 0; i < NCELLS; ++i)
        if (s->solid[i]) pix[i] = body;
    for (int k = 0; k < af->npts; ++k) {
        int k2 = (k + 1) % af->npts;
        draw_line_f(pix, af->px[k], af->py[k], af->px[k2], af->py[k2], edge);
    }
}

void render_arrows(const Sim *s, uint32_t *pix, int spacing)
{
    const uint32_t col = 0x00121216;
    for (int y = spacing / 2; y < NY; y += spacing) {
        for (int x = spacing / 2; x < NX; x += spacing) {
            int i = y * NX + x;
            if (s->solid[i]) continue;
            float vx = s->ux[i], vy = s->uy[i];
            float sp = sqrtf(vx * vx + vy * vy);
            if (sp < 0.004f) continue;

            float len = 24.f * sp / s->u0;
            if (len > 30.f) len = 30.f;
            float nx = vx / sp, ny = vy / sp;
            float x1 = (float)x + nx * len, y1 = (float)y + ny * len;
            draw_line_f(pix, (float)x, (float)y, x1, y1, col);

            /* Arrowhead: two barbs swept back +-30 degrees from the shaft. */
            float hx = -nx * 4.5f, hy = -ny * 4.5f;
            const float c30 = 0.866f, s30 = 0.5f;
            draw_line_f(pix, x1, y1, x1 + hx * c30 - hy * s30, y1 + hx * s30 + hy * c30, col);
            draw_line_f(pix, x1, y1, x1 + hx * c30 + hy * s30, y1 - hx * s30 + hy * c30, col);
        }
    }
}

/* tracer particles */

static unsigned g_rng = 0x12345u;
static float frand(void)
{
    /* xorshift32 -- fast, no libc state, deterministic */
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return (float)(g_rng & 0xFFFFFF) / 16777216.f;
}

void particles_init(Particle *p, int n, unsigned seed)
{
    g_rng = seed | 1u;
    for (int k = 0; k < n; ++k) {
        p[k].x = frand() * (NX * 0.95f);
        p[k].y = 1.f + frand() * (NY - 2);
    }
}

/* Bilinear velocity sample that ignores solid corners. */
static void vel_at(const Sim *s, float x, float y, float *vx, float *vy)
{
    if (x < 0.f) x = 0.f;
    if (x > NX - 1.001f) x = NX - 1.001f;
    if (y < 0.f) y = 0.f;
    if (y > NY - 1.001f) y = NY - 1.001f;
    int   x0 = (int)x, y0 = (int)y;
    float fx = x - x0, fy = y - y0;
    int   i00 = y0 * NX + x0, i10 = i00 + 1, i01 = i00 + NX, i11 = i01 + 1;

    float w00 = (1.f - fx) * (1.f - fy), w10 = fx * (1.f - fy);
    float w01 = (1.f - fx) * fy,         w11 = fx * fy;
    if (s->solid[i00]) w00 = 0.f;
    if (s->solid[i10]) w10 = 0.f;
    if (s->solid[i01]) w01 = 0.f;
    if (s->solid[i11]) w11 = 0.f;
    float ws = w00 + w10 + w01 + w11;
    if (ws < 1e-6f) { *vx = 0.f; *vy = 0.f; return; }
    float inv = 1.f / ws;
    *vx = (w00 * s->ux[i00] + w10 * s->ux[i10] + w01 * s->ux[i01] + w11 * s->ux[i11]) * inv;
    *vy = (w00 * s->uy[i00] + w10 * s->uy[i10] + w01 * s->uy[i01] + w11 * s->uy[i11]) * inv;
}

void particles_update(const Sim *s, Particle *p, int n, float dt)
{
    for (int k = 0; k < n; ++k) {
        float vx, vy;
        vel_at(s, p[k].x, p[k].y, &vx, &vy);
        /* midpoint (RK2) integration */
        float mx = p[k].x + 0.5f * dt * vx;
        float my = p[k].y + 0.5f * dt * vy;
        vel_at(s, mx, my, &vx, &vy);
        p[k].x += dt * vx;
        p[k].y += dt * vy;

        int xi = (int)p[k].x, yi = (int)p[k].y;
        int dead = p[k].x < 0.f || p[k].x >= NX - 2 || p[k].y < 1.f || p[k].y >= NY - 1;
        if (!dead && s->solid[yi * NX + xi]) dead = 1;
        if (dead) {
            p[k].x = 1.f + frand() * 6.f;
            p[k].y = 1.f + frand() * (NY - 2);
        }
    }
}

void render_particles(const Sim *s, const Particle *p, int n, uint32_t *pix)
{
    (void)s;
    for (int k = 0; k < n; ++k) {
        int x = (int)p[k].x, y = (int)p[k].y;
        if (x < 0 || x >= NX - 1 || y < 0 || y >= NY - 1) continue;
        for (int dy = 0; dy <= 1; ++dy) {
            for (int dx = 0; dx <= 1; ++dx) {
                uint32_t *q = &pix[(y + dy) * NX + (x + dx)];
                uint32_t  c = *q;
                /* 60% blend toward white */
                uint32_t r = ((c >> 16 & 0xFF) * 2 + 255 * 3) / 5;
                uint32_t g = ((c >> 8  & 0xFF) * 2 + 255 * 3) / 5;
                uint32_t b = ((c       & 0xFF) * 2 + 255 * 3) / 5;
                *q = (r << 16) | (g << 8) | b;
            }
        }
    }
}

void render_legend(uint32_t *pix, int mode)
{
    const int w = 220, h = 12, x0 = 14, y0 = 14; /* grid coords: bottom-left */
    for (int x = -1; x <= w; ++x) {
        for (int y = -1; y <= h; ++y) {
            int gx = x0 + x, gy = y0 + y;
            if (gx < 0 || gx >= NX || gy < 0 || gy >= NY) continue;
            if (x < 0 || x == w || y < 0 || y == h) {
                pix[gy * NX + gx] = 0x00202020; /* border */
            } else {
                float t = (float)x / (float)(w - 1);
                pix[gy * NX + gx] = (mode == VIS_VORTICITY) ? cmap_div(t) : cmap_jet(t);
            }
        }
    }
}
