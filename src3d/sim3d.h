/*
 * AeroSim3D -- 3D open-air flow around a finite wing.
 *
 * Solver: D3Q19 Lattice Boltzmann, BGK collision + Smagorinsky LES.
 * Axes: x streamwise (inflow at x=0 moving +x), y vertical, z spanwise.
 * Lattice units dx = dt = 1, cs^2 = 1/3. y/z boundaries are periodic
 * behind absorbing far-field layers, so the domain behaves as open air.
 */
#ifndef SIM3D_H
#define SIM3D_H

#include <stdint.h>
#include <stdbool.h>

#define NX3 224
#define NY3 96
#define NZ3 112
#define NCELLS3 (NX3 * NY3 * NZ3)
#define Q3 19

/* Absorbing far-field layer widths (cells) and peak blend strengths. */
#define ABC3_L  28
#define ABC3_R  52
#define ABC3_Y  22
#define ABC3_Z  20
#define ABC3_SMAX_SIDE 0.05f
#define ABC3_SMAX_OUT  0.10f

#define IDX3(x, y, z) (((z) * NY3 + (y)) * NX3 + (x))

typedef struct {
    float *f, *fnew;            /* SoA: f[q * NCELLS3 + cell] */
    uint8_t *solid, *careful;
    float *rho, *ux, *uy, *uz;  /* macroscopic fields          */
    float *sigma;               /* absorbing-layer strength    */

    float  u0;
    float  tau0;
    float  smag2;
    float  chord, span;         /* reference lengths, cells    */
    double refS;                /* reference area for CL/CD (cells^2) */
    double reynolds;

    long long step;
    int ramp_steps;

    double fx, fy, fz;          /* instantaneous force on the wing */
    double cl, cd;              /* smoothed coefficients (S = chord*span) */
} Sim3;

/* D3Q19 stencil tables, filled by sim3_init. */
extern int   E3[Q3][3];
extern int   OPP3[Q3];
extern float W3[Q3];

bool  sim3_init(Sim3 *s);
void  sim3_free(Sim3 *s);
void  sim3_reset_flow(Sim3 *s);
void  sim3_set_reynolds(Sim3 *s, double re);
void  sim3_step(Sim3 *s);
void  sim3_update_careful_mask(Sim3 *s);
void  sim3_refresh_cell(Sim3 *s, int i);
float sim3_inlet_speed(const Sim3 *s);

/* Q-criterion (vortex-core measure) of every interior cell into qout
 * (NCELLS3 floats): Q = 0.5 * (|Omega|^2 - |S|^2), > 0 inside vortices. */
void sim3_qcriterion(const Sim3 *s, float *qout);

#endif
