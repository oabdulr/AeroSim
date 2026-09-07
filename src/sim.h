/*
 * AeroSim -- 2D wind-tunnel simulation around an airfoil.
 *
 * Solver: D2Q9 Lattice Boltzmann Method (LBM) with BGK collision and a
 * Smagorinsky large-eddy-simulation (LES) subgrid model for stability and
 * physically plausible turbulence at high Reynolds numbers.
 *
 * Conventions:
 *   - Lattice units: dx = dt = 1, cs^2 = 1/3.
 *   - Grid is y-up (row 0 is the bottom of the tunnel). Flow enters at the
 *     left boundary moving in +x.
 *   - Kinematic viscosity nu = (tau - 1/2) / 3, Reynolds = u0 * chord / nu.
 */
#ifndef SIM_H
#define SIM_H

#include <stdint.h>
#include <stdbool.h>

/* Grid resolution in cells. Compile-time sized so the hot loop can use
 * constant strides. */
#define NX 960
#define NY 480
#define NCELLS (NX * NY)
#define QN 9

/* Width (in columns) of the viscous sponge zone that damps vortices before
 * they reach the outlet, suppressing pressure-wave reflections. */
#define SPONGE_W 50

/* Absorbing far-field layers (open-air boundary treatment). Inside these
 * layers the post-collision state is blended toward the undisturbed
 * freestream, so outgoing pressure waves and vortices are absorbed instead
 * of reflecting: the domain behaves like open air, not a closed tunnel. */
#define ABC_L 40   /* inlet-side layer width, cells   */
#define ABC_R 80   /* outlet-side layer width          */
#define ABC_TB 48  /* top and bottom layer width       */
#define ABC_SMAX_SIDE 0.05f
#define ABC_SMAX_OUT  0.10f

typedef struct {
    /* D2Q9 populations, double-buffered, SoA layout: f[q * NCELLS + cell]. */
    float *f, *fnew;

    /* Geometry masks. */
    uint8_t *solid;   /* 1 = cell inside the airfoil                     */
    uint8_t *careful; /* 1 = gather must check walls/boundaries/solids   */

    /* Macroscopic fields, refreshed every step (renderer + particles). */
    float *rho, *ux, *uy;

    /* Per-column relaxation time; includes the outlet sponge ramp. */
    float *tau_col;

    /* Absorbing-layer strength per cell (0 in the interior). */
    float *sigma;

    /* Parameters. */
    float  u0;      /* freestream speed, lattice units (kept ~0.1 for low Ma) */
    float  tau0;    /* base BGK relaxation time, derived from Reynolds        */
    float  smag2;   /* Smagorinsky constant squared, Cs^2                     */
    float  chord;   /* airfoil chord in cells (reference length for Re/Cl/Cd) */
    double reynolds;

    /* State. */
    long long step;
    int ramp_steps; /* inlet speed ramps 0 -> u0 over this many steps */

    /* Aerodynamic force on the airfoil via momentum exchange (lattice units). */
    double fx, fy;  /* instantaneous, this step        */
    double cl, cd;  /* exponentially smoothed coefficients */
} Sim;

bool  sim_init(Sim *s);
void  sim_free(Sim *s);
void  sim_reset_flow(Sim *s);
void  sim_set_reynolds(Sim *s, double re);
void  sim_step(Sim *s);
void  sim_update_careful_mask(Sim *s);
void  sim_refresh_cell(Sim *s, int i); /* set a cell to rest equilibrium */
float sim_inlet_speed(const Sim *s);   /* current (ramped) inlet speed   */

#endif
