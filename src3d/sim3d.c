#include "sim3d.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT restrict
#endif

/* D3Q19: rest, 6 axis directions, 12 edge diagonals. */
int   E3[Q3][3];
int   OPP3[Q3];
float W3[Q3];
static int OFF3[Q3]; /* linear memory offset of e_q */

static void build_tables(void)
{
    static const int e[Q3][3] = {
        { 0, 0, 0},
        { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
        { 1, 1, 0}, {-1,-1, 0}, { 1,-1, 0}, {-1, 1, 0},
        { 1, 0, 1}, {-1, 0,-1}, { 1, 0,-1}, {-1, 0, 1},
        { 0, 1, 1}, { 0,-1,-1}, { 0, 1,-1}, { 0,-1, 1},
    };
    for (int q = 0; q < Q3; ++q) {
        E3[q][0] = e[q][0]; E3[q][1] = e[q][1]; E3[q][2] = e[q][2];
        int n = abs(e[q][0]) + abs(e[q][1]) + abs(e[q][2]);
        W3[q] = n == 0 ? 1.f / 3.f : (n == 1 ? 1.f / 18.f : 1.f / 36.f);
        OFF3[q] = (e[q][2] * NY3 + e[q][1]) * NX3 + e[q][0];
    }
    for (int q = 0; q < Q3; ++q) {
        OPP3[q] = -1;
        for (int p = 0; p < Q3; ++p)
            if (e[p][0] == -e[q][0] && e[p][1] == -e[q][1] && e[p][2] == -e[q][2])
                { OPP3[q] = p; break; }
    }
}

static inline float feq3(int q, float rho, float ux, float uy, float uz)
{
    float eu = E3[q][0] * ux + E3[q][1] * uy + E3[q][2] * uz;
    float u2 = ux * ux + uy * uy + uz * uz;
    return W3[q] * rho * (1.f + 3.f * eu + 4.5f * eu * eu - 1.5f * u2);
}

bool sim3_init(Sim3 *s)
{
    build_tables();
    memset(s, 0, sizeof *s);
    s->f       = (float *)malloc(sizeof(float) * (size_t)Q3 * NCELLS3);
    s->fnew    = (float *)malloc(sizeof(float) * (size_t)Q3 * NCELLS3);
    s->solid   = (uint8_t *)calloc(NCELLS3, 1);
    s->careful = (uint8_t *)calloc(NCELLS3, 1);
    s->rho     = (float *)malloc(sizeof(float) * NCELLS3);
    s->ux      = (float *)malloc(sizeof(float) * NCELLS3);
    s->uy      = (float *)malloc(sizeof(float) * NCELLS3);
    s->uz      = (float *)malloc(sizeof(float) * NCELLS3);
    s->sigma   = (float *)malloc(sizeof(float) * NCELLS3);
    if (!s->f || !s->fnew || !s->solid || !s->careful ||
        !s->rho || !s->ux || !s->uy || !s->uz || !s->sigma) {
        sim3_free(s);
        return false;
    }

    /* Absorbing far-field layers on all six faces (quadratic ramps). */
    for (int z = 0; z < NZ3; ++z)
        for (int y = 0; y < NY3; ++y)
            for (int x = 0; x < NX3; ++x) {
                float sg = 0.f, r;
                if (x < ABC3_L) {
                    r = (float)(ABC3_L - x) / ABC3_L;
                    if (ABC3_SMAX_SIDE * r * r > sg) sg = ABC3_SMAX_SIDE * r * r;
                }
                if (x >= NX3 - ABC3_R) {
                    r = (float)(x - (NX3 - ABC3_R)) / ABC3_R;
                    if (ABC3_SMAX_OUT * r * r > sg) sg = ABC3_SMAX_OUT * r * r;
                }
                if (y < ABC3_Y)        { r = (float)(ABC3_Y - y) / ABC3_Y;         if (ABC3_SMAX_SIDE*r*r > sg) sg = ABC3_SMAX_SIDE*r*r; }
                if (y >= NY3 - ABC3_Y) { r = (float)(y - (NY3 - ABC3_Y)) / ABC3_Y; if (ABC3_SMAX_SIDE*r*r > sg) sg = ABC3_SMAX_SIDE*r*r; }
                if (z < ABC3_Z)        { r = (float)(ABC3_Z - z) / ABC3_Z;         if (ABC3_SMAX_SIDE*r*r > sg) sg = ABC3_SMAX_SIDE*r*r; }
                if (z >= NZ3 - ABC3_Z) { r = (float)(z - (NZ3 - ABC3_Z)) / ABC3_Z; if (ABC3_SMAX_SIDE*r*r > sg) sg = ABC3_SMAX_SIDE*r*r; }
                s->sigma[IDX3(x, y, z)] = sg;
            }

    s->u0         = 0.09f;
    s->smag2      = 0.14f * 0.14f;  /* Cs = 0.14 for the coarser 3D grid */
    s->chord      = 36.f;
    s->span       = 48.f;
    s->refS       = 36.0 * 48.0;
    s->ramp_steps = 1200;
    sim3_set_reynolds(s, 5000.0);
    sim3_update_careful_mask(s);
    sim3_reset_flow(s);
    return true;
}

void sim3_free(Sim3 *s)
{
    free(s->f);     free(s->fnew);
    free(s->solid); free(s->careful);
    free(s->rho);   free(s->ux); free(s->uy); free(s->uz);
    free(s->sigma);
    memset(s, 0, sizeof *s);
}

void sim3_set_reynolds(Sim3 *s, double re)
{
    s->reynolds = re;
    float nu = s->u0 * s->chord / (float)re;
    float t  = 3.f * nu + 0.5f;
    if (t < 0.5008f) t = 0.5008f;
    s->tau0 = t;
}

void sim3_refresh_cell(Sim3 *s, int i)
{
    for (int q = 0; q < Q3; ++q) {
        s->f[(size_t)q * NCELLS3 + i]    = W3[q];
        s->fnew[(size_t)q * NCELLS3 + i] = W3[q];
    }
    s->rho[i] = 1.f;
    s->ux[i] = s->uy[i] = s->uz[i] = 0.f;
}

void sim3_reset_flow(Sim3 *s)
{
    s->step = 0;
    s->cl = s->cd = 0.0;
    s->fx = s->fy = s->fz = 0.0;
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < NCELLS3; ++i) {
        s->rho[i] = 1.f;
        s->ux[i] = s->uy[i] = s->uz[i] = 0.f;
        for (int q = 0; q < Q3; ++q) {
            s->f[(size_t)q * NCELLS3 + i]    = W3[q];
            s->fnew[(size_t)q * NCELLS3 + i] = W3[q];
        }
    }
}

void sim3_update_careful_mask(Sim3 *s)
{
    int z;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (z = 0; z < NZ3; ++z) {
        for (int y = 0; y < NY3; ++y) {
            for (int x = 0; x < NX3; ++x) {
                int i = IDX3(x, y, z);
                if (x == 0 || x == NX3 - 1 || y == 0 || y == NY3 - 1 ||
                    z == 0 || z == NZ3 - 1) {
                    s->careful[i] = 1;
                    continue;
                }
                s->careful[i] = 0;
                if (s->solid[i]) continue;
                for (int q = 1; q < Q3; ++q) {
                    if (s->solid[i - OFF3[q]]) { s->careful[i] = 1; break; }
                }
            }
        }
    }
}

float sim3_inlet_speed(const Sim3 *s)
{
    if (s->step >= s->ramp_steps) return s->u0;
    float t = (float)s->step / (float)s->ramp_steps;
    return s->u0 * 0.5f * (1.f - cosf(3.14159265f * t));
}

void sim3_step(Sim3 *s)
{
    const float uin = sim3_inlet_speed(s);
    const float *RESTRICT f  = s->f;
    float       *RESTRICT fn = s->fnew;
    const uint8_t *RESTRICT solid   = s->solid;
    const uint8_t *RESTRICT careful = s->careful;
    const float   *RESTRICT sigma   = s->sigma;
    float *RESTRICT rho_o = s->rho;
    float *RESTRICT ux_o = s->ux, *RESTRICT uy_o = s->uy, *RESTRICT uz_o = s->uz;
    const float smag2 = s->smag2, tau0 = s->tau0;
    double fxa = 0.0, fya = 0.0, fza = 0.0;

    float fw[Q3];
    for (int q = 0; q < Q3; ++q) fw[q] = feq3(q, 1.f, uin, 0.f, 0.f);

    int z;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:fxa,fya,fza)
#endif
    for (z = 0; z < NZ3; ++z) {
        for (int y = 0; y < NY3; ++y) {
            for (int x = 0; x < NX3; ++x) {
                const int i = IDX3(x, y, z);
                if (solid[i]) continue;

                float fin[Q3];
                if (!careful[i]) {
                    for (int q = 0; q < Q3; ++q)
                        fin[q] = f[(size_t)q * NCELLS3 + i - OFF3[q]];
                } else {
                    for (int q = 0; q < Q3; ++q) {
                        int xs = x - E3[q][0];
                        int ys = y - E3[q][1];
                        int zs = z - E3[q][2];
                        if (xs < 0 || xs >= NX3) {
                            /* inlet/outlet planes are rewritten below */
                            fin[q] = f[(size_t)q * NCELLS3 + i];
                            continue;
                        }
                        /* periodic wrap on y/z (behind the absorbing layer) */
                        if (ys < 0) ys += NY3; else if (ys >= NY3) ys -= NY3;
                        if (zs < 0) zs += NZ3; else if (zs >= NZ3) zs -= NZ3;
                        int j = IDX3(xs, ys, zs);
                        if (solid[j]) {
                            /* half-way bounce-back + momentum exchange */
                            float fv = f[(size_t)OPP3[q] * NCELLS3 + i];
                            fin[q] = fv;
                            fxa -= 2.0 * (double)E3[q][0] * fv;
                            fya -= 2.0 * (double)E3[q][1] * fv;
                            fza -= 2.0 * (double)E3[q][2] * fv;
                        } else {
                            fin[q] = f[(size_t)q * NCELLS3 + j];
                        }
                    }
                }

                float rho = 0.f, mx = 0.f, my = 0.f, mz = 0.f;
                for (int q = 0; q < Q3; ++q) {
                    rho += fin[q];
                    mx  += E3[q][0] * fin[q];
                    my  += E3[q][1] * fin[q];
                    mz  += E3[q][2] * fin[q];
                }
                float inv = 1.f / rho;
                float ux = mx * inv, uy = my * inv, uz = mz * inv;

                float um2 = ux * ux + uy * uy + uz * uz;
                if (um2 > 0.15f) {
                    float sc = sqrtf(0.15f / um2);
                    ux *= sc; uy *= sc; uz *= sc;
                }

                float fe[Q3];
                float pxx = 0.f, pyy = 0.f, pzz = 0.f;
                float pxy = 0.f, pxz = 0.f, pyz = 0.f;
                for (int q = 0; q < Q3; ++q) {
                    fe[q] = feq3(q, rho, ux, uy, uz);
                    float fneq = fin[q] - fe[q];
                    pxx += E3[q][0] * E3[q][0] * fneq;
                    pyy += E3[q][1] * E3[q][1] * fneq;
                    pzz += E3[q][2] * E3[q][2] * fneq;
                    pxy += E3[q][0] * E3[q][1] * fneq;
                    pxz += E3[q][0] * E3[q][2] * fneq;
                    pyz += E3[q][1] * E3[q][2] * fneq;
                }
                float qn = sqrtf(pxx * pxx + pyy * pyy + pzz * pzz
                                 + 2.f * (pxy * pxy + pxz * pxz + pyz * pyz));
                float taue = 0.5f * (tau0 + sqrtf(tau0 * tau0 + 25.4558441f * smag2 * qn * inv));
                float om = 1.f / taue;

                const float sg = sigma[i];
                if (sg > 0.f) {
                    for (int q = 0; q < Q3; ++q) {
                        float v = fin[q] - om * (fin[q] - fe[q]);
                        fn[(size_t)q * NCELLS3 + i] = v - sg * (v - fw[q]);
                    }
                } else {
                    for (int q = 0; q < Q3; ++q)
                        fn[(size_t)q * NCELLS3 + i] = fin[q] - om * (fin[q] - fe[q]);
                }

                rho_o[i] = rho;
                ux_o[i] = ux; uy_o[i] = uy; uz_o[i] = uz;
            }
        }
    }

    /* Inlet plane x = 0: prescribed velocity, neighbour density. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (z = 0; z < NZ3; ++z) {
        for (int y = 0; y < NY3; ++y) {
            int i = IDX3(0, y, z);
            float r = rho_o[i + 1];
            if (!(r > 0.5f && r < 2.f)) r = 1.f;
            for (int q = 0; q < Q3; ++q)
                fn[(size_t)q * NCELLS3 + i] = feq3(q, r, uin, 0.f, 0.f);
            rho_o[i] = r; ux_o[i] = uin; uy_o[i] = 0.f; uz_o[i] = 0.f;
        }
    }
    /* Outlet plane x = NX3-1: zeroth-order extrapolation. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (z = 0; z < NZ3; ++z) {
        for (int y = 0; y < NY3; ++y) {
            int i = IDX3(NX3 - 1, y, z);
            for (int q = 0; q < Q3; ++q)
                fn[(size_t)q * NCELLS3 + i] = fn[(size_t)q * NCELLS3 + i - 1];
            rho_o[i] = rho_o[i - 1];
            ux_o[i] = ux_o[i - 1]; uy_o[i] = uy_o[i - 1]; uz_o[i] = uz_o[i - 1];
        }
    }

    float *tmp = s->f; s->f = s->fnew; s->fnew = tmp;

    s->fx = fxa; s->fy = fya; s->fz = fza;
    float  uref = uin > 0.02f ? uin : 0.02f;
    double qd   = 0.5 * (double)uref * uref * s->refS;
    double a    = 0.004;
    s->cl += a * (fya / qd - s->cl);
    s->cd += a * (fxa / qd - s->cd);

    s->step++;
}

void sim3_qcriterion(const Sim3 *s, float *qout)
{
    const float *ux = s->ux, *uy = s->uy, *uz = s->uz;
    int z;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (z = 0; z < NZ3; ++z) {
        for (int y = 0; y < NY3; ++y) {
            for (int x = 0; x < NX3; ++x) {
                int i = IDX3(x, y, z);
                if (x < 1 || x >= NX3 - 1 || y < 1 || y >= NY3 - 1 ||
                    z < 1 || z >= NZ3 - 1 || s->solid[i]) {
                    qout[i] = 0.f;
                    continue;
                }
                const int sx = 1, sy = NX3, sz = NX3 * NY3;
                float dux = 0.5f * (ux[i+sx] - ux[i-sx]);
                float duy = 0.5f * (ux[i+sy] - ux[i-sy]);
                float duz = 0.5f * (ux[i+sz] - ux[i-sz]);
                float dvx = 0.5f * (uy[i+sx] - uy[i-sx]);
                float dvy = 0.5f * (uy[i+sy] - uy[i-sy]);
                float dvz = 0.5f * (uy[i+sz] - uy[i-sz]);
                float dwx = 0.5f * (uz[i+sx] - uz[i-sx]);
                float dwy = 0.5f * (uz[i+sy] - uz[i-sy]);
                float dwz = 0.5f * (uz[i+sz] - uz[i-sz]);
                /* Q = 0.5(|Omega|^2 - |S|^2); positive marks vortex cores. */
                float sxx = dux, syy = dvy, szz = dwz;
                float sxy = 0.5f * (duy + dvx);
                float sxz = 0.5f * (duz + dwx);
                float syz = 0.5f * (dvz + dwy);
                float oxy = 0.5f * (duy - dvx);
                float oxz = 0.5f * (duz - dwx);
                float oyz = 0.5f * (dvz - dwy);
                float S2 = sxx*sxx + syy*syy + szz*szz + 2.f*(sxy*sxy + sxz*sxz + syz*syz);
                float O2 = 2.f * (oxy*oxy + oxz*oxz + oyz*oyz);
                qout[i] = 0.5f * (O2 - S2);
            }
        }
    }
}
