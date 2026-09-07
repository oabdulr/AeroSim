#include "sim.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT restrict
#endif

/* D2Q9 stencil.
 *      6  2  5
 *      3  0  1        EX/EY: lattice velocities
 *      7  4  8        OPP:   opposite direction (for bounce-back)
 *                     MIRY:  mirrored about the x axis (for free-slip walls)
 */
static const int   EX[QN]   = { 0, 1, 0,-1, 0, 1,-1,-1, 1 };
static const int   EY[QN]   = { 0, 0, 1, 0,-1, 1, 1,-1,-1 };
static const int   OPP[QN]  = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
static const int   MIRY[QN] = { 0, 1, 4, 3, 2, 8, 7, 6, 5 };
static const float W[QN] = { 4.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f,
                             1.f/36.f, 1.f/36.f, 1.f/36.f, 1.f/36.f };

/* Second-order Maxwell-Boltzmann equilibrium. */
static inline float feq(int q, float rho, float ux, float uy)
{
    float eu = (float)EX[q] * ux + (float)EY[q] * uy;
    float u2 = ux * ux + uy * uy;
    return W[q] * rho * (1.f + 3.f * eu + 4.5f * eu * eu - 1.5f * u2);
}

bool sim_init(Sim *s)
{
    memset(s, 0, sizeof *s);
    s->f       = (float *)malloc(sizeof(float) * (size_t)QN * NCELLS);
    s->fnew    = (float *)malloc(sizeof(float) * (size_t)QN * NCELLS);
    s->solid   = (uint8_t *)calloc(NCELLS, 1);
    s->careful = (uint8_t *)calloc(NCELLS, 1);
    s->rho     = (float *)malloc(sizeof(float) * NCELLS);
    s->ux      = (float *)malloc(sizeof(float) * NCELLS);
    s->uy      = (float *)malloc(sizeof(float) * NCELLS);
    s->tau_col = (float *)malloc(sizeof(float) * NX);
    s->sigma   = (float *)malloc(sizeof(float) * NCELLS);
    if (!s->f || !s->fnew || !s->solid || !s->careful ||
        !s->rho || !s->ux || !s->uy || !s->tau_col || !s->sigma) {
        sim_free(s);
        return false;
    }

    /* Absorbing far-field layers: quadratic ramp of blend strength toward
     * each edge; a cell in a corner takes the strongest of the overlapping
     * ramps. Zero in the interior, so the hot loop is unaffected there. */
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            float sg = 0.f, r;
            if (x < ABC_L) {
                r = (float)(ABC_L - x) / ABC_L;
                if (ABC_SMAX_SIDE * r * r > sg) sg = ABC_SMAX_SIDE * r * r;
            }
            if (x >= NX - ABC_R) {
                r = (float)(x - (NX - ABC_R)) / ABC_R;
                if (ABC_SMAX_OUT * r * r > sg) sg = ABC_SMAX_OUT * r * r;
            }
            if (y < ABC_TB) {
                r = (float)(ABC_TB - y) / ABC_TB;
                if (ABC_SMAX_SIDE * r * r > sg) sg = ABC_SMAX_SIDE * r * r;
            }
            if (y >= NY - ABC_TB) {
                r = (float)(y - (NY - ABC_TB)) / ABC_TB;
                if (ABC_SMAX_SIDE * r * r > sg) sg = ABC_SMAX_SIDE * r * r;
            }
            s->sigma[y * NX + x] = sg;
        }
    }

    s->u0         = 0.10f;
    s->smag2      = 0.13f * 0.13f;
    s->chord      = 150.f;
    s->ramp_steps = 2500;
    sim_set_reynolds(s, 10000.0);
    sim_update_careful_mask(s);
    sim_reset_flow(s);
    return true;
}

void sim_free(Sim *s)
{
    free(s->f);      free(s->fnew);
    free(s->solid);  free(s->careful);
    free(s->rho);    free(s->ux);     free(s->uy);
    free(s->tau_col); free(s->sigma);
    memset(s, 0, sizeof *s);
}

void sim_set_reynolds(Sim *s, double re)
{
    s->reynolds = re;
    float nu = s->u0 * s->chord / (float)re;
    float t  = 3.f * nu + 0.5f;
    if (t < 0.5005f) t = 0.5005f; /* hard stability floor */
    s->tau0 = t;

    /* Outlet sponge: viscosity rises quadratically over the last columns. */
    for (int x = 0; x < NX; ++x) {
        float tt = t;
        int   x0 = NX - SPONGE_W;
        if (x >= x0) {
            float sf = (float)(x - x0) / (float)SPONGE_W;
            tt = t + (1.6f - t) * sf * sf;
        }
        s->tau_col[x] = tt;
    }
}

void sim_refresh_cell(Sim *s, int i)
{
    for (int q = 0; q < QN; ++q) {
        s->f[(size_t)q * NCELLS + i]    = W[q];
        s->fnew[(size_t)q * NCELLS + i] = W[q];
    }
    s->rho[i] = 1.f;
    s->ux[i]  = 0.f;
    s->uy[i]  = 0.f;
}

void sim_reset_flow(Sim *s)
{
    s->step = 0;
    s->cl = s->cd = 0.0;
    s->fx = s->fy = 0.0;
    for (int i = 0; i < NCELLS; ++i) {
        s->rho[i] = 1.f;
        s->ux[i]  = 0.f;
        s->uy[i]  = 0.f;
    }
    for (int q = 0; q < QN; ++q) {
        float *a = s->f    + (size_t)q * NCELLS;
        float *b = s->fnew + (size_t)q * NCELLS;
        for (int i = 0; i < NCELLS; ++i) { a[i] = W[q]; b[i] = W[q]; }
    }
}

void sim_update_careful_mask(Sim *s)
{
    memset(s->careful, 0, NCELLS);
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = y * NX + x;
            if (x == 0 || x == NX - 1 || y == 0 || y == NY - 1) {
                s->careful[i] = 1;
                continue;
            }
            if (s->solid[i]) continue;
            for (int q = 1; q < QN; ++q) {
                if (s->solid[(y - EY[q]) * NX + (x - EX[q])]) {
                    s->careful[i] = 1;
                    break;
                }
            }
        }
    }
}

float sim_inlet_speed(const Sim *s)
{
    if (s->step >= s->ramp_steps) return s->u0;
    float t = (float)s->step / (float)s->ramp_steps;
    return s->u0 * 0.5f * (1.f - cosf(3.14159265f * t));
}

/* One LBM time step: fused pull-streaming + collision over the whole grid,
 * then inlet/outlet columns, then buffer swap and force bookkeeping. */
void sim_step(Sim *s)
{
    const float           uin     = sim_inlet_speed(s);
    const float *RESTRICT f       = s->f;
    float       *RESTRICT fn      = s->fnew;
    const uint8_t *RESTRICT solid   = s->solid;
    const uint8_t *RESTRICT careful = s->careful;
    const float  *RESTRICT tau_col  = s->tau_col;
    const float  *RESTRICT sigma    = s->sigma;
    float *RESTRICT rho_o = s->rho;
    float *RESTRICT ux_o  = s->ux;
    float *RESTRICT uy_o  = s->uy;
    const float smag2 = s->smag2;
    double fxa = 0.0, fya = 0.0;

    /* Freestream equilibrium at the current (ramped) inlet speed -- the
     * target state the absorbing layers relax toward. */
    float fw[QN];
    for (int q = 0; q < QN; ++q) fw[q] = feq(q, 1.f, uin, 0.f);
    int y; /* MSVC's OpenMP wants the parallel loop index declared out here */

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:fxa,fya)
#endif
    for (y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            const int i = y * NX + x;
            if (solid[i]) continue;

            float fin[QN];
            if (!careful[i]) {
                /* Fast path: all eight neighbours are interior fluid. */
                fin[0] = f[0 * (size_t)NCELLS + i];
                fin[1] = f[1 * (size_t)NCELLS + i - 1];
                fin[2] = f[2 * (size_t)NCELLS + i - NX];
                fin[3] = f[3 * (size_t)NCELLS + i + 1];
                fin[4] = f[4 * (size_t)NCELLS + i + NX];
                fin[5] = f[5 * (size_t)NCELLS + i - NX - 1];
                fin[6] = f[6 * (size_t)NCELLS + i - NX + 1];
                fin[7] = f[7 * (size_t)NCELLS + i + NX + 1];
                fin[8] = f[8 * (size_t)NCELLS + i + NX - 1];
            } else {
                for (int q = 0; q < QN; ++q) {
                    int xs = x - EX[q], ys = y - EY[q];
                    if (ys < 0 || ys >= NY) {
                        /* Free-slip tunnel wall: specular reflection. The
                         * mirrored population leaves from this row. */
                        int xc = xs < 0 ? 0 : (xs >= NX ? NX - 1 : xs);
                        fin[q] = f[(size_t)MIRY[q] * NCELLS + (size_t)y * NX + xc];
                    } else if (xs < 0 || xs >= NX) {
                        /* Inlet/outlet columns are rewritten below. */
                        fin[q] = f[(size_t)q * NCELLS + i];
                    } else {
                        int j = ys * NX + xs;
                        if (solid[j]) {
                            /* Half-way bounce-back off the airfoil, plus
                             * momentum exchange: the reflected population was
                             * travelling with e_opp = -e_q into the wall. */
                            float fv = f[(size_t)OPP[q] * NCELLS + i];
                            fin[q] = fv;
                            fxa -= 2.0 * (double)EX[q] * fv;
                            fya -= 2.0 * (double)EY[q] * fv;
                        } else {
                            fin[q] = f[(size_t)q * NCELLS + j];
                        }
                    }
                }
            }

            /* Macroscopic moments. */
            float rho = fin[0] + fin[1] + fin[2] + fin[3] + fin[4]
                      + fin[5] + fin[6] + fin[7] + fin[8];
            float inv = 1.f / rho;
            float ux = (fin[1] - fin[3] + fin[5] - fin[6] - fin[7] + fin[8]) * inv;
            float uy = (fin[2] - fin[4] + fin[5] + fin[6] - fin[7] - fin[8]) * inv;

            /* Safety valve: cap |u| well below lattice sound speed so one
             * pathological cell cannot detonate the whole field. */
            float um2 = ux * ux + uy * uy;
            if (um2 > 0.16f) {
                float sc = sqrtf(0.16f / um2);
                ux *= sc; uy *= sc;
            }

            /* Equilibria and the non-equilibrium stress tensor Pi. */
            float fe[QN];
            float pxx = 0.f, pyy = 0.f, pxy = 0.f;
            for (int q = 0; q < QN; ++q) {
                fe[q] = feq(q, rho, ux, uy);
                float fneq = fin[q] - fe[q];
                pxx += (float)(EX[q] * EX[q]) * fneq;
                pyy += (float)(EY[q] * EY[q]) * fneq;
                pxy += (float)(EX[q] * EY[q]) * fneq;
            }

            /* Smagorinsky LES: closed-form effective relaxation time
             * tau_eff = (tau + sqrt(tau^2 + 18*sqrt(2)*Cs^2*|Pi|/rho)) / 2. */
            float qn   = sqrtf(pxx * pxx + pyy * pyy + 2.f * pxy * pxy);
            float tau  = tau_col[x];
            float taue = 0.5f * (tau + sqrtf(tau * tau + 25.4558441f * smag2 * qn * inv));
            float om   = 1.f / taue;

            const float sg = sigma[i];
            if (sg > 0.f) {
                /* Absorbing far-field layer: relax toward the freestream. */
                for (int q = 0; q < QN; ++q) {
                    float v = fin[q] - om * (fin[q] - fe[q]);
                    fn[(size_t)q * NCELLS + i] = v - sg * (v - fw[q]);
                }
            } else {
                for (int q = 0; q < QN; ++q)
                    fn[(size_t)q * NCELLS + i] = fin[q] - om * (fin[q] - fe[q]);
            }

            rho_o[i] = rho;
            ux_o[i]  = ux;
            uy_o[i]  = uy;
        }
    }

    /* Inlet (x = 0): prescribed velocity, density taken from the neighbour
     * column so pressure waves can pass through. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (y = 0; y < NY; ++y) {
        int   i = y * NX;
        float r = rho_o[i + 1];
        if (!(r > 0.5f && r < 2.f)) r = 1.f;
        for (int q = 0; q < QN; ++q)
            fn[(size_t)q * NCELLS + i] = feq(q, r, uin, 0.f);
        rho_o[i] = r; ux_o[i] = uin; uy_o[i] = 0.f;
    }

    /* Outlet (x = NX-1): zeroth-order extrapolation behind the sponge. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (y = 0; y < NY; ++y) {
        int i = y * NX + NX - 1;
        for (int q = 0; q < QN; ++q)
            fn[(size_t)q * NCELLS + i] = fn[(size_t)q * NCELLS + i - 1];
        rho_o[i] = rho_o[i - 1]; ux_o[i] = ux_o[i - 1]; uy_o[i] = uy_o[i - 1];
    }

    /* Swap buffers. */
    float *tmp = s->f; s->f = s->fnew; s->fnew = tmp;

    /* Force -> smoothed coefficients. Reference dynamic pressure uses the
     * current inlet speed (floored during early ramp to avoid divide-by-0). */
    s->fx = fxa;
    s->fy = fya;
    float  uref = uin > 0.02f ? uin : 0.02f;
    double qd   = 0.5 * (double)uref * uref * s->chord;
    double a    = 0.002;
    s->cl += a * (fya / qd - s->cl);
    s->cd += a * (fxa / qd - s->cd);

    s->step++;
}
