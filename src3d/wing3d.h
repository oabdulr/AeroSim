/* Finite-span NACA wing with a part-span plain flap. */
#ifndef WING3D_H
#define WING3D_H

#include "sim3d.h"

#define WG_MAX_PTS 256

typedef struct {
    float camber, camber_pos, thickness; /* NACA 4-digit, fractions of chord */
    float aoa_deg;                       /* + = nose up                      */
    float flap_deg;                      /* + = trailing edge down           */
    float flap_hinge;                    /* hinge x/c                        */
    float flap_frac;                     /* flapped fraction of span (inboard) */
    float chord, span;                   /* cells                            */
    float cx, cy, cz;                    /* quarter-chord mid-span position  */
} Wing3;

/* Section outline at a given spanwise station (grid x-y coords, y-up).
 * flapped != 0 applies the flap deflection. Returns point count. */
int wing3_section(const Wing3 *w, int flapped, float *px, float *py);

/* Rasterize the wing into the sim's solid mask (and refresh changed cells). */
void wing3_rasterize(const Wing3 *w, Sim3 *s);

/* Triangle mesh for rendering: interleaved [nx ny nz  x y z] per vertex,
 * 3 vertices per triangle. Returns the vertex count written (<= maxverts). */
int wing3_mesh(const Wing3 *w, float *out, int maxverts);

/* Detailed render mesh: separated slotted flap with cove and gap,
 * flap-track fairings, hydraulic actuators that extend with deflection,
 * vortex generators, panel lines baked into matcol (one RGB per vertex).
 * Purely visual -- the solid mask stays the sealed aerodynamic wing. */
int wing3_mesh_detailed(const Wing3 *w, float *out, float *matcol, int maxverts);

/* Light anchor points: 0 = port nav (red), 1 = starboard nav (green),
 * 2/3 = tip strobes (white). Fills pos[4][3] in grid coordinates. */
void wing3_lights(const Wing3 *w, float pos[4][3]);

#endif
