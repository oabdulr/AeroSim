/* AeroSim3D -- interactive 3D LBM flow around a finite wing.
 * Win32 + OpenGL 1.1 (both ship with Windows), no external dependencies. */
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <windowsx.h>
#include <commdlg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sim3d.h"
#include "wing3d.h"
#include "gl3d.h"
#include "gpu3d.h"
#include "model3d.h"

/* state */
static Sim3   g_sim;
static Wing3  g_wing;
static Camera3 g_cam;
static HWND   g_hwnd;
static int    g_cw = 1280, g_ch = 800;
static BOOL   g_run = TRUE, g_paused = FALSE, g_help = TRUE;
static int    g_spf = 2;

static int    g_mode = 0;          /* slice field: 0 Cp, 1 |u|, 2 vorticity */
static int    g_orient = 0;        /* slice plane: 0 z, 1 y, 2 x            */
static int    g_slicepos[3];       /* per-orientation slice coordinate      */
static BOOL   g_slice = FALSE, g_splats = TRUE, g_particles = FALSE;
static BOOL   g_stream = TRUE, g_hairs = TRUE;
static BOOL   g_gpu = FALSE;
static float  g_qscale = 24.f;     /* Q-criterion threshold factor          */
static float  g_flap_mem = 20.f;
static double g_fps, g_mlups;

/* wing/model render mesh (+ per-vertex surface pressure colors) */
#define MESH_MAX 660000
static float *g_mesh;
static float *g_meshcol;   /* colors actually drawn                    */
static float *g_meshmat;   /* baked material/livery colors             */
static int    g_meshn;
static BOOL   g_paint = FALSE; /* FALSE = Cp colormap, TRUE = livery   */

/* navigation lights: port red, starboard green, two white strobes */
static float  g_lpos[4][3];
static double g_time;

/* aero benchmark table: press K to record the current body's numbers */
#define BENCH_MAX 6
typedef struct { char name[36]; float aoa; double cd, cl; } Bench;
static Bench g_bench[BENCH_MAX];
static int   g_benchn;

/* surface flow hairs (2 verts per hair) */
#define HAIR_MAX 4000
static float *g_hpos, *g_hcol;
static int    g_hn;

/* streaklines: smoke trails emitted from a fixed rake. Each trail point is
 * a particle advected by the flow, so a passing vortex bends only the part
 * of the trail inside it -- unlike instantaneous streamlines, which re-route
 * wholesale every frame. */
#define SK_N   63          /* 9 spanwise x 7 vertical seeds */
#define SK_LEN 300         /* points per trail: ~0.7 cells of travel per
                              frame, so 300 spans the whole test section */
static float   *g_skp;     /* positions, SK_N * SK_LEN * 3  */
static float   *g_skspd;   /* local speed per point         */
static uint8_t *g_skalive;
static float   *g_slpos, *g_slcol; /* per-run draw scratch  */

/* Q-criterion splats */
#define SPLAT_MAX 400000
static float *g_qval, *g_spos, *g_scol;
static int    g_sn;

/* tracer particles */
#define NP3 15000
typedef struct { float x, y, z; } P3;
static P3 g_p3[NP3];
static float *g_ppos, *g_pcol;

/* slice texture: 256 x 128 RGBA */
#define STW 256
#define STH 128
static uint8_t g_stex[STW * STH * 4];

/* small helpers */
static unsigned g_rng = 77777u;
static float frand3(void)
{
    g_rng ^= g_rng << 13; g_rng ^= g_rng >> 17; g_rng ^= g_rng << 5;
    return (float)(g_rng & 0xFFFFFF) / 16777216.f;
}

static inline float clamp01(float t) { return t < 0.f ? 0.f : (t > 1.f ? 1.f : t); }

static void jetf(float t, float *r, float *g, float *b)
{
    t = clamp01(t);
    *r = clamp01(1.5f - fabsf(4.f * t - 3.f));
    *g = clamp01(1.5f - fabsf(4.f * t - 2.f));
    *b = clamp01(1.5f - fabsf(4.f * t - 1.f));
    if (t < 0.125f) *b = 0.55f + 3.6f * t;
    if (t > 0.875f) *r = 1.f - (t - 0.875f) * 1.6f;
}

static void divf(float t, float *r, float *g, float *b)
{
    t = clamp01(t);
    if (t < 0.5f) {
        float u = t * 2.f;
        *r = 0.12f + 0.85f * u; *g = 0.22f + 0.75f * u; *b = 0.83f + 0.14f * u;
    } else {
        float u = (t - 0.5f) * 2.f;
        *r = 0.97f - 0.18f * u; *g = 0.97f - 0.88f * u; *b = 0.97f - 0.85f * u;
    }
}

static const char *mode_name(int m)
{
    return m == 0 ? "Pressure (Cp)" : (m == 1 ? "Velocity |u|" : "Vorticity");
}
static const char *orient_name(int o)
{
    return o == 0 ? "spanwise (side view)" : (o == 1 ? "horizontal (top view)" : "crossflow (front view)");
}

static uint8_t *g_oldsolid, *g_resetmask;
static float g_myaw; /* model yaw, degrees */

static void rebuild_wing(void)
{
    if (g_oldsolid) memcpy(g_oldsolid, g_sim.solid, NCELLS3);
    if (model3_active()) {
        model3_rasterize(&g_sim, g_wing.aoa_deg, g_myaw);
        g_meshn = model3_mesh(g_mesh, MESH_MAX, g_wing.aoa_deg, g_myaw);
        g_sim.refS = model3_refS();
        for (int v = 0; v < g_meshn; ++v) { /* neutral livery for models */
            g_meshmat[v*3+0] = 0.80f;
            g_meshmat[v*3+1] = 0.81f;
            g_meshmat[v*3+2] = 0.84f;
        }
    } else {
        wing3_rasterize(&g_wing, &g_sim);
        g_meshn = wing3_mesh_detailed(&g_wing, g_mesh, g_meshmat, MESH_MAX);
        g_sim.refS = (double)g_wing.chord * g_wing.span;
        wing3_lights(&g_wing, g_lpos);
    }
    if (g_gpu && g_oldsolid && g_resetmask) {
        for (int i = 0; i < NCELLS3; ++i)
            g_resetmask[i] = (uint8_t)(g_oldsolid[i] && !g_sim.solid[i]);
        gpu3_update_geometry(&g_sim, g_resetmask);
    }
}

/* slice sampling */
static float g_rhoref = 1.f;

static void update_rhoref(void)
{
    /* undisturbed upstream strip inside the absorbing layer */
    double m = 0.0; int n = 0;
    for (int z = NZ3 / 2 - 16; z < NZ3 / 2 + 16; ++z)
        for (int y = NY3 / 2 - 16; y < NY3 / 2 + 16; ++y)
            for (int x = 4; x < 10; ++x) { m += g_sim.rho[IDX3(x, y, z)]; ++n; }
    g_rhoref = (float)(m / n);
}

static void field_color(int i, int di_a, int di_b, uint8_t *px)
{
    /* di_a/di_b: strides of the two in-plane axes, for the vorticity curl */
    float r, g, b;
    if (g_sim.solid[i]) { px[0] = 235; px[1] = 235; px[2] = 240; px[3] = 255; return; }
    if (g_mode == 0) {
        float cp = ((g_sim.rho[i] - g_rhoref) / 3.f) / (0.5f * g_sim.u0 * g_sim.u0);
        jetf((cp + 2.f) / 3.f, &r, &g, &b);
    } else if (g_mode == 1) {
        float vx = g_sim.ux[i], vy = g_sim.uy[i], vz = g_sim.uz[i];
        jetf(sqrtf(vx * vx + vy * vy + vz * vz) / (1.6f * g_sim.u0), &r, &g, &b);
    } else {
        /* vorticity component perpendicular to the slice plane */
        float w;
        if (g_orient == 0)      /* z-normal: dv/dx - du/dy  */
            w = 0.5f * (g_sim.uy[i + 1] - g_sim.uy[i - 1])
              - 0.5f * (g_sim.ux[i + NX3] - g_sim.ux[i - NX3]);
        else if (g_orient == 1) /* y-normal: du/dz - dw/dx  */
            w = 0.5f * (g_sim.ux[i + NX3 * NY3] - g_sim.ux[i - NX3 * NY3])
              - 0.5f * (g_sim.uz[i + 1] - g_sim.uz[i - 1]);
        else                    /* x-normal: dw/dy - dv/dz  */
            w = 0.5f * (g_sim.uz[i + NX3] - g_sim.uz[i - NX3])
              - 0.5f * (g_sim.uy[i + NX3 * NY3] - g_sim.uy[i - NX3 * NY3]);
        (void)di_a; (void)di_b;
        float ws = 18.f * g_sim.u0 / g_sim.chord;
        divf(0.5f + 0.5f * w / ws, &r, &g, &b);
    }
    px[0] = (uint8_t)(r * 255.f); px[1] = (uint8_t)(g * 255.f);
    px[2] = (uint8_t)(b * 255.f); px[3] = 255;
}

static void build_slice(float *umax, float *vmax,
                        float c0[3], float c1[3], float c2[3], float c3[3])
{
    memset(g_stex, 0, sizeof g_stex);
    if (g_mode == 0) update_rhoref();

    if (g_orient == 0) {
        int zs = g_slicepos[0];
        for (int y = 0; y < NY3; ++y)
            for (int x = 0; x < NX3; ++x) {
                int i = IDX3(x, y, zs);
                if (x < 1 || x >= NX3-1 || y < 1 || y >= NY3-1) continue;
                field_color(i, 1, NX3, &g_stex[(y * STW + x) * 4]);
            }
        *umax = (float)NX3 / STW; *vmax = (float)NY3 / STH;
        float z = (float)zs;
        c0[0]=0;   c0[1]=0;   c0[2]=z;
        c1[0]=NX3; c1[1]=0;   c1[2]=z;
        c2[0]=NX3; c2[1]=NY3; c2[2]=z;
        c3[0]=0;   c3[1]=NY3; c3[2]=z;
    } else if (g_orient == 1) {
        int ys = g_slicepos[1];
        for (int z = 0; z < NZ3; ++z)
            for (int x = 0; x < NX3; ++x) {
                int i = IDX3(x, ys, z);
                if (x < 1 || x >= NX3-1 || z < 1 || z >= NZ3-1) continue;
                field_color(i, 1, NX3 * NY3, &g_stex[(z * STW + x) * 4]);
            }
        *umax = (float)NX3 / STW; *vmax = (float)NZ3 / STH;
        float y = (float)ys;
        c0[0]=0;   c0[1]=y; c0[2]=0;
        c1[0]=NX3; c1[1]=y; c1[2]=0;
        c2[0]=NX3; c2[1]=y; c2[2]=NZ3;
        c3[0]=0;   c3[1]=y; c3[2]=NZ3;
    } else {
        int xs = g_slicepos[2];
        for (int y = 0; y < NY3; ++y)
            for (int z = 0; z < NZ3; ++z) {
                int i = IDX3(xs, y, z);
                if (z < 1 || z >= NZ3-1 || y < 1 || y >= NY3-1) continue;
                field_color(i, NX3, NX3 * NY3, &g_stex[(y * STW + z) * 4]);
            }
        *umax = (float)NZ3 / STW; *vmax = (float)NY3 / STH;
        float x = (float)xs;
        c0[0]=x; c0[1]=0;   c0[2]=0;
        c1[0]=x; c1[1]=0;   c1[2]=NZ3;
        c2[0]=x; c2[1]=NY3; c2[2]=NZ3;
        c3[0]=x; c3[1]=NY3; c3[2]=0;
    }
}

/* Q-criterion splats */
static void build_splats(void)
{
    sim3_qcriterion(&g_sim, g_qval);
    float uc = g_sim.u0 / g_sim.chord;
    float thr = g_qscale * uc * uc;
    g_sn = 0;
    for (int z = 2; z < NZ3 - 2; z += 1) {
        for (int y = 2; y < NY3 - 2; y += 1) {
            const float *qrow = &g_qval[IDX3(0, y, z)];
            for (int x = 2; x < NX3 - 2; x += 1) {
                if (qrow[x] <= thr) continue;
                if (g_sn >= SPLAT_MAX) return;
                int i = IDX3(x, y, z);
                float vx = g_sim.ux[i], vy = g_sim.uy[i], vz = g_sim.uz[i];
                float sp = sqrtf(vx * vx + vy * vy + vz * vz) / (1.5f * g_sim.u0);
                float r, g, b;
                jetf(sp, &r, &g, &b);
                /* stronger vortices glow brighter */
                float a = 0.06f + 0.20f * clamp01(qrow[x] / (8.f * thr));
                g_spos[g_sn * 3 + 0] = (float)x;
                g_spos[g_sn * 3 + 1] = (float)y;
                g_spos[g_sn * 3 + 2] = (float)z;
                g_scol[g_sn * 4 + 0] = r;
                g_scol[g_sn * 4 + 1] = g;
                g_scol[g_sn * 4 + 2] = b;
                g_scol[g_sn * 4 + 3] = a;
                ++g_sn;
            }
        }
    }
}

/* particles */
static void p3_respawn(P3 *p)
{
    p->x = 2.f + frand3() * 8.f;
    p->y = ABC3_Y + frand3() * (NY3 - 2.f * ABC3_Y);
    p->z = ABC3_Z + frand3() * (NZ3 - 2.f * ABC3_Z);
}

static void vel3_at(float x, float y, float z, float *vx, float *vy, float *vz)
{
    if (x < 0) x = 0; if (x > NX3 - 1.001f) x = NX3 - 1.001f;
    if (y < 0) y = 0; if (y > NY3 - 1.001f) y = NY3 - 1.001f;
    if (z < 0) z = 0; if (z > NZ3 - 1.001f) z = NZ3 - 1.001f;
    int x0 = (int)x, y0 = (int)y, z0 = (int)z;
    float fx = x - x0, fy = y - y0, fz = z - z0;
    float ax = 0, ay = 0, az = 0, ws = 0;
    for (int dz = 0; dz <= 1; ++dz)
        for (int dy = 0; dy <= 1; ++dy)
            for (int dx = 0; dx <= 1; ++dx) {
                int i = IDX3(x0 + dx, y0 + dy, z0 + dz);
                if (g_sim.solid[i]) continue;
                float w = (dx ? fx : 1.f - fx) * (dy ? fy : 1.f - fy) * (dz ? fz : 1.f - fz);
                ax += w * g_sim.ux[i]; ay += w * g_sim.uy[i]; az += w * g_sim.uz[i];
                ws += w;
            }
    if (ws < 1e-6f) { *vx = *vy = *vz = 0.f; return; }
    *vx = ax / ws; *vy = ay / ws; *vz = az / ws;
}

static void particles3_update(float dt)
{
    for (int k = 0; k < NP3; ++k) {
        P3 *p = &g_p3[k];
        float vx, vy, vz;
        vel3_at(p->x, p->y, p->z, &vx, &vy, &vz);
        vel3_at(p->x + 0.5f * dt * vx, p->y + 0.5f * dt * vy, p->z + 0.5f * dt * vz,
                &vx, &vy, &vz);
        p->x += dt * vx; p->y += dt * vy; p->z += dt * vz;
        int bad = p->x < 0 || p->x >= NX3 - 2 || p->y < 1 || p->y >= NY3 - 1 ||
                  p->z < 1 || p->z >= NZ3 - 1;
        if (!bad && g_sim.solid[IDX3((int)p->x, (int)p->y, (int)p->z)]) bad = 1;
        if (bad) p3_respawn(p);
        g_ppos[k * 3 + 0] = p->x;
        g_ppos[k * 3 + 1] = p->y;
        g_ppos[k * 3 + 2] = p->z;
    }
}

/* surface pressure coloring + flow hairs (cow-style) */
static float rho3_at(float x, float y, float z)
{
    if (x < 0) x = 0; if (x > NX3 - 1.001f) x = NX3 - 1.001f;
    if (y < 0) y = 0; if (y > NY3 - 1.001f) y = NY3 - 1.001f;
    if (z < 0) z = 0; if (z > NZ3 - 1.001f) z = NZ3 - 1.001f;
    int x0 = (int)x, y0 = (int)y, z0 = (int)z;
    float fx = x - x0, fy = y - y0, fz = z - z0;
    float acc = 0.f, ws = 0.f;
    for (int dz = 0; dz <= 1; ++dz)
        for (int dy = 0; dy <= 1; ++dy)
            for (int dx = 0; dx <= 1; ++dx) {
                int i = IDX3(x0 + dx, y0 + dy, z0 + dz);
                if (g_sim.solid[i]) continue;
                float w = (dx ? fx : 1.f - fx) * (dy ? fy : 1.f - fy) * (dz ? fz : 1.f - fz);
                acc += w * g_sim.rho[i]; ws += w;
            }
    return ws > 1e-6f ? acc / ws : g_rhoref;
}

static void build_surface_colors(void)
{
    if (g_paint) { /* livery mode: baked material colors */
        memcpy(g_meshcol, g_meshmat, sizeof(float) * 3 * (size_t)g_meshn);
        return;
    }
    /* one field sample per triangle (imported models can have 200k+ tris) */
    const float inv_qd = 1.f / (0.5f * g_sim.u0 * g_sim.u0);
    const int ntri = g_meshn / 3;
    int t;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (t = 0; t < ntri; ++t) {
        const float *n = &g_mesh[t * 18];
        const float *p = n + 3;
        float rho = rho3_at(p[0] + n[0] * 2.2f, p[1] + n[1] * 2.2f, p[2] + n[2] * 2.2f);
        float cp = ((rho - g_rhoref) / 3.f) * inv_qd;
        float r, g, b;
        jetf((cp + 2.f) / 3.f, &r, &g, &b);
        for (int v = 0; v < 3; ++v) {
            g_meshcol[(t * 3 + v) * 3 + 0] = r;
            g_meshcol[(t * 3 + v) * 3 + 1] = g;
            g_meshcol[(t * 3 + v) * 3 + 2] = b;
        }
    }
}

static void build_hairs(void)
{
    g_hn = 0;
    int stride = g_meshn / 3600;
    if (stride < 21) stride = 21;
    for (int v = 0; v < g_meshn && g_hn < HAIR_MAX; v += stride) {
        const float *n = &g_mesh[v * 6];
        const float *p = n + 3;
        float bx = p[0] + n[0] * 1.6f, by = p[1] + n[1] * 1.6f, bz = p[2] + n[2] * 1.6f;
        float vx, vy, vz;
        vel3_at(bx + n[0] * 1.2f, by + n[1] * 1.2f, bz + n[2] * 1.2f, &vx, &vy, &vz);
        float dot = vx * n[0] + vy * n[1] + vz * n[2];
        float tx = vx - dot * n[0], ty = vy - dot * n[1], tz = vz - dot * n[2];
        float sp = sqrtf(tx * tx + ty * ty + tz * tz);
        if (sp < 0.15f * g_sim.u0) continue;
        float len = 3.f + 9.f * clamp01(sp / g_sim.u0);
        float inv = len / sp;
        int k = g_hn * 6;
        g_hpos[k + 0] = bx;             g_hpos[k + 1] = by;             g_hpos[k + 2] = bz;
        g_hpos[k + 3] = bx + tx * inv;  g_hpos[k + 4] = by + ty * inv;  g_hpos[k + 5] = bz + tz * inv;
        int c = g_hn * 8;
        for (int e = 0; e < 2; ++e) {
            g_hcol[c + e * 4 + 0] = 0.05f;
            g_hcol[c + e * 4 + 1] = 0.05f;
            g_hcol[c + e * 4 + 2] = 0.09f;
            g_hcol[c + e * 4 + 3] = e == 0 ? 0.9f : 0.55f;
        }
        ++g_hn;
    }
}

/* streaklines (smoke rake) */
static void seed_at(int s, float *x, float *y, float *z)
{
    int zi = s / 7 - 4;
    int yi = s % 7 - 3;
    *x = (float)(ABC3_L + 2); /* just clear of the inlet absorbing layer */
    *y = g_wing.cy + yi * 8.5f;
    *z = g_wing.cz + zi * (g_wing.span * 0.145f);
}

static void streaks_clear(void)
{
    if (g_skalive) memset(g_skalive, 0, (size_t)SK_N * SK_LEN);
}

static void streaks_update(float dt)
{
    /* Emit a fresh trail point only after the freestream has advanced a
     * fixed distance, so trail length and spacing are independent of the
     * steps-per-frame setting (at 1 step/frame trails used to collapse
     * into invisible stubs). */
    static float acc = 0.f;
    acc += g_sim.u0 * dt;
    int emit = acc >= 0.75f;
    if (emit) acc = 0.f;

    for (int s = 0; s < SK_N; ++s) {
        float   *P  = &g_skp[(size_t)s * SK_LEN * 3];
        float   *SP = &g_skspd[(size_t)s * SK_LEN];
        uint8_t *A  = &g_skalive[(size_t)s * SK_LEN];

        /* advect every live trail point (midpoint RK2) */
        for (int k = 0; k < SK_LEN; ++k) {
            if (!A[k]) continue;
            float x = P[k*3], y = P[k*3+1], z = P[k*3+2];
            float vx, vy, vz;
            vel3_at(x, y, z, &vx, &vy, &vz);
            vel3_at(x + 0.5f*dt*vx, y + 0.5f*dt*vy, z + 0.5f*dt*vz, &vx, &vy, &vz);
            x += dt * vx; y += dt * vy; z += dt * vz;
            if (x < 1.f || x >= NX3 - 2.f || y < 1.f || y >= NY3 - 1.f ||
                z < 1.f || z >= NZ3 - 1.f ||
                g_sim.solid[IDX3((int)x, (int)y, (int)z)]) {
                A[k] = 0;
                continue;
            }
            P[k*3] = x; P[k*3+1] = y; P[k*3+2] = z;
            SP[k] = sqrtf(vx*vx + vy*vy + vz*vz);
        }

        /* age the trail and emit a fresh point at the seed */
        if (emit) {
            memmove(&P[3],  &P[0],  sizeof(float) * 3 * (SK_LEN - 1));
            memmove(&SP[1], &SP[0], sizeof(float) * (SK_LEN - 1));
            memmove(&A[1],  &A[0],  SK_LEN - 1);
            seed_at(s, &P[0], &P[1], &P[2]);
            float vx, vy, vz;
            vel3_at(P[0], P[1], P[2], &vx, &vy, &vz);
            SP[0] = sqrtf(vx*vx + vy*vy + vz*vz);
            A[0] = 1;
        }
    }
}

static void streaks_draw(void)
{
    /* Real smoke dissipates: opacity decays with age, and segments that
     * turbulence has stretched hard get dimmer and finally tear. Without
     * this, old folded filaments pile up into visual spaghetti. */
    const float u0 = g_sim.u0;
    for (int s = 0; s < SK_N; ++s) {
        const float   *P  = &g_skp[(size_t)s * SK_LEN * 3];
        const float   *SP = &g_skspd[(size_t)s * SK_LEN];
        const uint8_t *A  = &g_skalive[(size_t)s * SK_LEN];
        int k = 0;
        while (k < SK_LEN) {
            if (!A[k]) { ++k; continue; }
            int m = 0;
            while (k < SK_LEN && A[k] && m < SK_LEN) {
                float x = P[k*3], y = P[k*3+1], z = P[k*3+2];
                float stretch = 1.f;
                if (m > 0) {
                    float dx = x - g_slpos[(m-1)*3];
                    float dy = y - g_slpos[(m-1)*3+1];
                    float dz = z - g_slpos[(m-1)*3+2];
                    float d = sqrtf(dx*dx + dy*dy + dz*dz);
                    if (d > 7.f) break;          /* filament torn apart */
                    if (d > 1.5f) stretch = 1.5f / d;
                }
                float age = expf((float)k * -(1.f / 95.f));
                float a = 0.85f * age * stretch;
                if (a < 0.025f) break;           /* fully dissipated */
                g_slpos[m*3+0] = x; g_slpos[m*3+1] = y; g_slpos[m*3+2] = z;
                float r, g, b;
                jetf(SP[k] / (1.5f * u0), &r, &g, &b);
                g_slcol[m*4+0] = r; g_slcol[m*4+1] = g; g_slcol[m*4+2] = b;
                g_slcol[m*4+3] = a;
                ++m; ++k;
            }
            if (m >= 2) gl3_draw_strip(g_slpos, g_slcol, m, 2.f);
            if (m == 0) ++k; /* ensure progress if the first point tore */
        }
    }
}

/* window */
static int g_dragging = 0, g_lastmx, g_lastmy;

static LRESULT CALLBACK wndproc3(HWND h, UINT m, WPARAM w, LPARAM l)
{
    switch (m) {
    case WM_SIZE:
        g_cw = LOWORD(l); g_ch = HIWORD(l);
        gl3_resize(g_cw, g_ch);
        return 0;
    case WM_ERASEBKGND: return 1;
    case WM_PAINT: { PAINTSTRUCT ps; BeginPaint(h, &ps); EndPaint(h, &ps); return 0; }
    case WM_CLOSE: case WM_DESTROY: g_run = FALSE; PostQuitMessage(0); return 0;
    case WM_LBUTTONDOWN:
        g_dragging = 1; g_lastmx = GET_X_LPARAM(l); g_lastmy = GET_Y_LPARAM(l);
        SetCapture(h);
        return 0;
    case WM_LBUTTONUP:
        g_dragging = 0; ReleaseCapture();
        return 0;
    case WM_MOUSEMOVE:
        if (g_dragging) {
            int mx = GET_X_LPARAM(l), my = GET_Y_LPARAM(l);
            g_cam.yaw   += 0.4f * (mx - g_lastmx);
            g_cam.pitch += 0.3f * (my - g_lastmy);
            if (g_cam.pitch > 89.f)  g_cam.pitch = 89.f;
            if (g_cam.pitch < -89.f) g_cam.pitch = -89.f;
            g_lastmx = mx; g_lastmy = my;
        }
        return 0;
    case WM_MOUSEWHEEL: {
        int d = GET_WHEEL_DELTA_WPARAM(w);
        g_cam.dist *= d > 0 ? 0.9f : 1.111f;
        if (g_cam.dist < 60.f)  g_cam.dist = 60.f;
        if (g_cam.dist > 900.f) g_cam.dist = 900.f;
        return 0;
    }
    case WM_KEYDOWN:
        switch (w) {
        case VK_ESCAPE: g_run = FALSE; PostQuitMessage(0); break;
        case VK_UP:
            g_wing.aoa_deg += 1.f; if (g_wing.aoa_deg > 25.f) g_wing.aoa_deg = 25.f;
            rebuild_wing(); break;
        case VK_DOWN:
            g_wing.aoa_deg -= 1.f; if (g_wing.aoa_deg < -25.f) g_wing.aoa_deg = -25.f;
            rebuild_wing(); break;
        case VK_LEFT: /* yaw imported models about the vertical axis */
            if (model3_active()) {
                g_myaw -= 5.f; if (g_myaw < -180.f) g_myaw += 360.f;
                rebuild_wing();
            }
            break;
        case VK_RIGHT:
            if (model3_active()) {
                g_myaw += 5.f; if (g_myaw > 180.f) g_myaw -= 360.f;
                rebuild_wing();
            }
            break;
        case 'F':
            if (model3_active()) break; /* the flap belongs to the wing */
            if (g_wing.flap_deg != 0.f) { g_flap_mem = g_wing.flap_deg; g_wing.flap_deg = 0.f; }
            else g_wing.flap_deg = g_flap_mem;
            rebuild_wing(); break;
        case VK_OEM_4: /* [ */
            if (model3_active()) break;
            g_wing.flap_deg -= 5.f; if (g_wing.flap_deg < -20.f) g_wing.flap_deg = -20.f;
            rebuild_wing(); break;
        case VK_OEM_6: /* ] */
            if (model3_active()) break;
            g_wing.flap_deg += 5.f; if (g_wing.flap_deg > 40.f) g_wing.flap_deg = 40.f;
            rebuild_wing(); break;
        case '1': g_mode = 0; g_slice = TRUE; break;
        case '2': g_mode = 1; g_slice = TRUE; break;
        case '3': g_mode = 2; g_slice = TRUE; break;
        case 'M': { /* import an STL/OBJ model */
            char fname[MAX_PATH]; fname[0] = 0;
            OPENFILENAMEA ofn;
            memset(&ofn, 0, sizeof ofn);
            ofn.lStructSize = sizeof ofn;
            ofn.hwndOwner   = h;
            ofn.lpstrFilter = "3D models (*.stl;*.obj)\0*.stl;*.obj\0All files\0*.*\0";
            ofn.lpstrFile   = fname;
            ofn.nMaxFile    = MAX_PATH;
            ofn.lpstrTitle  = "Import 3D model into the flow";
            ofn.Flags       = OFN_FILEMUSTEXIST | OFN_NOCHANGEDIR;
            if (GetOpenFileNameA(&ofn)) {
                char err[256];
                int ok = 0;
#if defined(_MSC_VER)
                /* belt-and-braces: a malformed file must never take the
                 * whole app down */
                __try {
                    ok = model3_load(fname, err, sizeof err);
                    if (ok) {
                        g_wing.aoa_deg = 0.f;
                        g_myaw = 0.f;
                        rebuild_wing();
                    }
                } __except (EXCEPTION_EXECUTE_HANDLER) {
                    ok = 0;
                    model3_clear();
                    snprintf(err, sizeof err,
                             "The file triggered a fault while importing "
                             "(is it a valid, closed mesh?)");
                    rebuild_wing();
                }
#else
                ok = model3_load(fname, err, sizeof err);
                if (ok) { g_wing.aoa_deg = 0.f; rebuild_wing(); }
#endif
                if (ok) {
                    sim3_reset_flow(&g_sim);
                    if (g_gpu) gpu3_reset_all(&g_sim);
                    streaks_clear();
                } else {
                    MessageBoxA(h, err, "Model import failed", MB_ICONWARNING);
                }
            }
            break;
        }
        case 'N': /* back to the built-in wing */
            if (model3_active()) {
                model3_clear();
                g_wing.aoa_deg = 8.f;
                rebuild_wing();
                sim3_reset_flow(&g_sim);
                if (g_gpu) gpu3_reset_all(&g_sim);
                streaks_clear();
            }
            break;
        case 'U': /* cycle model axis convention (z-up files etc.) */
            if (model3_active()) {
                model3_cycle_axes();
                rebuild_wing();
                sim3_reset_flow(&g_sim);
                if (g_gpu) gpu3_reset_all(&g_sim);
                streaks_clear();
            }
            break;
        case 'O': g_orient = (g_orient + 1) % 3; break;
        case 'Z': {
            int lim = g_orient == 0 ? NZ3 : (g_orient == 1 ? NY3 : NX3);
            g_slicepos[g_orient] -= 2;
            if (g_slicepos[g_orient] < 2) g_slicepos[g_orient] = 2;
            (void)lim; break;
        }
        case 'X': {
            int lim = g_orient == 0 ? NZ3 : (g_orient == 1 ? NY3 : NX3);
            g_slicepos[g_orient] += 2;
            if (g_slicepos[g_orient] > lim - 3) g_slicepos[g_orient] = lim - 3;
            break;
        }
        case 'L': g_slice = !g_slice; break;
        case 'V': g_splats = !g_splats; break;
        case 'P': g_particles = !g_particles; break;
        case 'S': g_stream = !g_stream; break;
        case 'T': g_hairs = !g_hairs; break;
        case 'C': g_paint = !g_paint; break;
        case 'H': g_help = !g_help; break;
        case 'R':
            sim3_reset_flow(&g_sim);
            if (g_gpu) gpu3_reset_all(&g_sim);
            streaks_clear();
            break;
        case VK_SPACE: g_paused = !g_paused; break;
        case '9': g_qscale *= 1.5f; if (g_qscale > 400.f) g_qscale = 400.f; break;
        case '0': g_qscale /= 1.5f; if (g_qscale < 0.5f)  g_qscale = 0.5f;  break;
        case '7': /* wind speed down */
            g_sim.u0 *= 0.90f;
            if (g_sim.u0 < 0.04f) g_sim.u0 = 0.04f;
            sim3_set_reynolds(&g_sim, g_sim.reynolds); /* retune viscosity */
            break;
        case '8': /* wind speed up */
            g_sim.u0 *= 1.10f;
            if (g_sim.u0 > 0.13f) g_sim.u0 = 0.13f;
            sim3_set_reynolds(&g_sim, g_sim.reynolds);
            break;
        case 'K': { /* record a benchmark entry for A/B comparisons */
            int slot = g_benchn % BENCH_MAX;
            snprintf(g_bench[slot].name, sizeof g_bench[slot].name, "%.28s%s",
                     model3_active() ? model3_name() : "NACA2412 wing",
                     model3_active() ? "" : (g_wing.flap_deg != 0.f ? "+flap" : ""));
            g_bench[slot].aoa = g_wing.aoa_deg;
            g_bench[slot].cd  = g_sim.cd;
            g_bench[slot].cl  = g_sim.cl;
            ++g_benchn;
            break;
        }
        case VK_OEM_COMMA: {
            double re = g_sim.reynolds / 1.5; if (re < 300) re = 300;
            sim3_set_reynolds(&g_sim, re); break;
        }
        case VK_OEM_PERIOD: {
            double re = g_sim.reynolds * 1.5; if (re > 30000) re = 30000;
            sim3_set_reynolds(&g_sim, re); break;
        }
        case VK_OEM_MINUS: case VK_SUBTRACT:
            g_spf -= 1; if (g_spf < 1) g_spf = 1; break;
        case VK_OEM_PLUS: case VK_ADD:
            g_spf += 1; if (g_spf > 24) g_spf = 24; break;
        }
        return 0;
    }
    return DefWindowProcA(h, m, w, l);
}

int WINAPI WinMain(HINSTANCE hi, HINSTANCE hp, LPSTR cmd, int show)
{
    (void)hp; (void)show;
    SetProcessDPIAware();

    WNDCLASSA wc;
    memset(&wc, 0, sizeof wc);
    wc.style = CS_OWNDC;
    wc.lpfnWndProc = wndproc3;
    wc.hInstance = hi;
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.lpszClassName = "AeroSim3DWnd";
    RegisterClassA(&wc);

    RECT wr = {0, 0, g_cw, g_ch};
    AdjustWindowRect(&wr, WS_OVERLAPPEDWINDOW, FALSE);
    g_hwnd = CreateWindowExA(0, "AeroSim3DWnd",
                             "AeroSim3D - finite wing, open air (D3Q19 LBM + LES)",
                             WS_OVERLAPPEDWINDOW | WS_VISIBLE,
                             CW_USEDEFAULT, CW_USEDEFAULT,
                             wr.right - wr.left, wr.bottom - wr.top,
                             NULL, NULL, hi, NULL);
    if (!g_hwnd) return 1;
    if (!gl3_init(g_hwnd)) {
        MessageBoxA(NULL, "OpenGL initialization failed.", "AeroSim3D", MB_ICONERROR);
        return 1;
    }
    gl3_resize(g_cw, g_ch);

    if (!sim3_init(&g_sim)) {
        MessageBoxA(NULL, "Out of memory allocating the 3D grid (~450 MB needed).",
                    "AeroSim3D", MB_ICONERROR);
        return 1;
    }
    g_mesh    = (float *)malloc(sizeof(float) * 6 * MESH_MAX);
    g_meshcol = (float *)malloc(sizeof(float) * 3 * MESH_MAX);
    g_meshmat = (float *)malloc(sizeof(float) * 3 * MESH_MAX);
    if (!g_meshmat) return 1;
    g_qval = (float *)malloc(sizeof(float) * NCELLS3);
    g_spos = (float *)malloc(sizeof(float) * 3 * SPLAT_MAX);
    g_scol = (float *)malloc(sizeof(float) * 4 * SPLAT_MAX);
    g_ppos = (float *)malloc(sizeof(float) * 3 * NP3);
    g_pcol = (float *)malloc(sizeof(float) * 4 * NP3);
    g_hpos = (float *)malloc(sizeof(float) * 6 * HAIR_MAX);
    g_hcol = (float *)malloc(sizeof(float) * 8 * HAIR_MAX);
    g_skp     = (float *)malloc(sizeof(float) * 3 * SK_N * SK_LEN);
    g_skspd   = (float *)malloc(sizeof(float) * SK_N * SK_LEN);
    g_skalive = (uint8_t *)calloc((size_t)SK_N * SK_LEN, 1);
    g_slpos = (float *)malloc(sizeof(float) * 3 * SK_LEN);
    g_slcol = (float *)malloc(sizeof(float) * 4 * SK_LEN);
    g_oldsolid  = (uint8_t *)malloc(NCELLS3);
    g_resetmask = (uint8_t *)malloc(NCELLS3);
    if (!g_mesh || !g_meshcol || !g_qval || !g_spos || !g_scol || !g_ppos ||
        !g_pcol || !g_hpos || !g_hcol || !g_skp || !g_skspd || !g_skalive ||
        !g_slpos || !g_slcol || !g_oldsolid || !g_resetmask) return 1;

    memset(&g_wing, 0, sizeof g_wing);
    g_wing.camber = 0.02f; g_wing.camber_pos = 0.4f; g_wing.thickness = 0.12f;
    g_wing.aoa_deg = 8.f;  g_wing.flap_deg = 0.f;
    g_wing.flap_hinge = 0.70f; g_wing.flap_frac = 0.60f;
    g_wing.chord = g_sim.chord; g_wing.span = g_sim.span;
    g_wing.cx = NX3 * 0.32f; g_wing.cy = NY3 * 0.52f; g_wing.cz = NZ3 * 0.5f;

    /* optional command-line model path */
    if (cmd && *cmd) {
        char path[MAX_PATH];
        snprintf(path, sizeof path, "%s", cmd);
        char *pp = path;
        if (*pp == '"') { ++pp; char *q = strrchr(pp, '"'); if (q) *q = 0; }
        char err[256];
        if (model3_load(pp, err, sizeof err)) { g_wing.aoa_deg = 0.f; g_myaw = 0.f; }
        else MessageBoxA(NULL, err, "Model load failed", MB_ICONWARNING);
    }
    rebuild_wing();

    /* GPU compute backend; silently fall back to the OpenMP CPU solver */
    g_gpu = gpu3_init(&g_sim);
    g_spf = g_gpu ? 8 : 2;

    g_slicepos[0] = (int)g_wing.cz;
    g_slicepos[1] = (int)g_wing.cy;
    g_slicepos[2] = (int)(g_wing.cx + g_wing.chord);

    g_cam.yaw = -38.f; g_cam.pitch = 18.f; g_cam.dist = 300.f;
    g_cam.tx = g_wing.cx + 20.f; g_cam.ty = g_wing.cy; g_cam.tz = g_wing.cz;

    for (int k = 0; k < NP3; ++k) {
        p3_respawn(&g_p3[k]);
        g_p3[k].x = frand3() * (NX3 * 0.9f); /* scatter through the domain */
        g_ppos[k * 3 + 0] = g_p3[k].x;
        g_ppos[k * 3 + 1] = g_p3[k].y;
        g_ppos[k * 3 + 2] = g_p3[k].z;
        g_pcol[k * 4 + 0] = 1.f; g_pcol[k * 4 + 1] = 1.f;
        g_pcol[k * 4 + 2] = 1.f; g_pcol[k * 4 + 3] = 0.45f;
    }

    LARGE_INTEGER freq, t0, t1;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&t0);

    MSG msg;
    while (g_run) {
        while (PeekMessageA(&msg, NULL, 0, 0, PM_REMOVE)) {
            if (msg.message == WM_QUIT) g_run = FALSE;
            TranslateMessage(&msg);
            DispatchMessageA(&msg);
        }
        if (!g_run) break;

        int steps = 0;
        if (!g_paused) {
            if (g_gpu) {
                gpu3_steps(&g_sim, g_spf);
                gpu3_readback(&g_sim, g_spf);
            } else {
                for (int k = 0; k < g_spf; ++k) sim3_step(&g_sim);
            }
            steps = g_spf;
            if (g_particles) particles3_update((float)steps);
        }

        /* fast-math-proof watchdog */
        float probe = g_sim.rho[IDX3(NX3 / 2, NY3 / 2, NZ3 / 2)];
        uint32_t pbits; memcpy(&pbits, &probe, sizeof pbits);
        if (((pbits >> 23) & 0xFFu) == 0xFFu || probe < 0.1f || probe > 10.f) {
            sim3_reset_flow(&g_sim);
            if (g_gpu) gpu3_reset_all(&g_sim);
        }

        update_rhoref();
        build_surface_colors();
        if (g_hairs) build_hairs(); else g_hn = 0;
        if (g_stream && steps > 0) streaks_update((float)steps);
        if (g_splats) build_splats(); else g_sn = 0;
        if (g_particles) {
            for (int k = 0; k < NP3; ++k) {
                float vx, vy, vz;
                vel3_at(g_ppos[k*3], g_ppos[k*3+1], g_ppos[k*3+2], &vx, &vy, &vz);
                float r, gg, b;
                jetf(sqrtf(vx*vx + vy*vy + vz*vz) / (1.5f * g_sim.u0), &r, &gg, &b);
                g_pcol[k*4+0] = r; g_pcol[k*4+1] = gg; g_pcol[k*4+2] = b;
                g_pcol[k*4+3] = 0.5f;
            }
        }

        gl3_begin_frame(&g_cam);
        gl3_draw_box(0, 0, 0, NX3, NY3, NZ3);
        gl3_draw_mesh_colored(g_mesh, g_meshcol, g_meshn);
        if (g_hn) gl3_draw_lines(g_hpos, g_hcol, g_hn, 1.f);
        if (g_stream) streaks_draw();
        if (g_slice) {
            float um, vm, c0[3], c1[3], c2[3], c3[3];
            build_slice(&um, &vm, c0, c1, c2, c3);
            gl3_draw_slice(g_stex, STW, STH, um, vm, c0, c1, c2, c3);
        }
        if (g_sn) gl3_draw_points(g_spos, g_scol, g_sn, 3.f);
        if (g_particles) gl3_draw_points(g_ppos, g_pcol, NP3, 1.5f);

        /* navigation lights + double-flash strobes (wing only) */
        if (!model3_active()) {
            float t = (float)fmod(g_time, 1.4);
            int strobe = (t < 0.06f) || (t > 0.18f && t < 0.24f);
            float sa = strobe ? 1.f : 0.f;
            float lc[4 * 4] = {
                1.f, 0.05f, 0.05f, 0.95f,   /* port red        */
                0.05f, 1.f, 0.10f, 0.95f,   /* starboard green */
                1.f, 1.f, 1.f, sa,
                1.f, 1.f, 1.f, sa,
            };
            float halo[4 * 4];
            for (int k = 0; k < 4; ++k) {
                halo[k*4+0] = lc[k*4+0]; halo[k*4+1] = lc[k*4+1];
                halo[k*4+2] = lc[k*4+2]; halo[k*4+3] = lc[k*4+3] * 0.20f;
            }
            gl3_draw_points(&g_lpos[0][0], lc, 4, 6.f);
            gl3_draw_points(&g_lpos[0][0], halo, 4, 17.f);
        }

        gl3_overlay_begin(g_cw, g_ch);
        char buf[320];
        int ty = 22;
        if (model3_active()) {
            snprintf(buf, sizeof buf,
                     "Model: %s (%d tris)   Pitch %+.0f   Yaw %+.0f (LEFT/RIGHT)   Re %.0f   wind %d%%   t* %.1f   |   %s %.3f   |   N = wing, U = fix axes",
                     model3_name(), model3_tris(), (double)g_wing.aoa_deg,
                     (double)g_myaw,
                     g_sim.reynolds, (int)(g_sim.u0 / 0.09f * 100.f + 0.5f),
                     (double)g_sim.step * g_sim.u0 / g_sim.chord,
                     g_sim.cl < -0.02 ? "DOWNFORCE coeff" : "lift coeff",
                     g_sim.cl < 0 ? -g_sim.cl : g_sim.cl);
        } else {
            snprintf(buf, sizeof buf,
                     "NACA 2412 wing  AoA %+.0f deg   Flap %+.0f deg (inner %d%% span)   Re %.0f   wind %d%%   t* %.1f",
                     (double)g_wing.aoa_deg, (double)g_wing.flap_deg,
                     (int)(g_wing.flap_frac * 100.f), g_sim.reynolds,
                     (int)(g_sim.u0 / 0.09f * 100.f + 0.5f),
                     (double)g_sim.step * g_sim.u0 / g_sim.chord);
        }
        gl3_overlay_text(12, ty, buf); ty += 19;
        snprintf(buf, sizeof buf,
                 "CL %+.3f   CD %+.4f   L/D %+.1f   |   engine: %s   |   %d steps/frame   %.0f fps   %.0f MLUPS%s%s",
                 g_sim.cl, g_sim.cd, fabs(g_sim.cd) > 1e-9 ? g_sim.cl / g_sim.cd : 0.0,
                 g_gpu ? "GPU compute (GL 4.3)" : "CPU (OpenMP)",
                 g_spf, g_fps, g_mlups,
                 g_paused ? "   [PAUSED]" : "",
                 g_sim.step < g_sim.ramp_steps ? "   [spinning up...]" : "");
        gl3_overlay_text(12, ty, buf); ty += 19;
        if (g_help) {
            gl3_overlay_text(12, ty,
                "drag = orbit   wheel = zoom   UP/DOWN AoA/pitch   F flap   [ ] flap angle   M import model (STL/OBJ)   C paint (Cp / livery)   S smoke   T hairs   V vortex cores   9/0 threshold");
            ty += 19;
            snprintf(buf, sizeof buf,
                "L slice (%s, %s; 1/2/3 field, O plane, Z/X move)   P particles   , . Reynolds   7/8 wind speed   K record benchmark   +/- sim speed   SPACE pause   R reset   H help",
                mode_name(g_mode), orient_name(g_orient));
            gl3_overlay_text(12, ty, buf); ty += 19;
        }
        if (g_benchn > 0) {
            ty += 6;
            gl3_overlay_text(12, ty, "--- benchmarks (K to record; lower CD = less drag, negative CL = downforce) ---");
            ty += 19;
            int nshow = g_benchn < BENCH_MAX ? g_benchn : BENCH_MAX;
            for (int k = 0; k < nshow; ++k) {
                snprintf(buf, sizeof buf, "%-30s  AoA %+5.1f   CD %6.4f   CL %+7.4f   L/D %+6.2f",
                         g_bench[k].name, (double)g_bench[k].aoa,
                         g_bench[k].cd, g_bench[k].cl,
                         fabs(g_bench[k].cd) > 1e-9 ? g_bench[k].cl / g_bench[k].cd : 0.0);
                gl3_overlay_text(12, ty, buf);
                ty += 17;
            }
        }
        gl3_overlay_end();
        gl3_end_frame();

        QueryPerformanceCounter(&t1);
        double dt = (double)(t1.QuadPart - t0.QuadPart) / (double)freq.QuadPart;
        t0 = t1;
        if (dt > 1e-6) {
            g_time  += dt;
            g_fps   += 0.1 * (1.0 / dt - g_fps);
            g_mlups += 0.1 * ((double)steps * NCELLS3 / dt / 1e6 - g_mlups);
        }
        if (g_paused) Sleep(10);
    }

    gl3_shutdown();
    sim3_free(&g_sim);
    return 0;
}
