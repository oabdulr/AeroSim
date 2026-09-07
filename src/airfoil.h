/* Parametric NACA 4-digit airfoil with a hinged plain flap, plus
 * rasterization of its outline into the simulation's solid mask. */
#ifndef AIRFOIL_H
#define AIRFOIL_H

#include "sim.h"

#define AF_MAX_PTS 512

typedef struct {
    /* NACA 4-digit shape parameters, as fractions of chord.
     * NACA 2412 -> camber 0.02, camber_pos 0.4, thickness 0.12. */
    float camber;
    float camber_pos;
    float thickness;

    float aoa_deg;    /* angle of attack, + = nose up            */
    float flap_deg;   /* flap deflection, + = trailing edge down */
    float flap_hinge; /* hinge chordwise position, x/c           */

    float chord;      /* chord length in grid cells */
    float cx, cy;     /* quarter-chord location in the grid (y-up) */

    /* Generated closed outline polygon in grid coordinates. */
    float px[AF_MAX_PTS], py[AF_MAX_PTS];
    int   npts;
} Airfoil;

void airfoil_build(Airfoil *af);                    /* recompute outline    */
void airfoil_rasterize(const Airfoil *af, Sim *s);  /* update solid mask    */

#endif
