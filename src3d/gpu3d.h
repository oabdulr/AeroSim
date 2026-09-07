/* GPU compute backend for the D3Q19 solver: the full stream-collide step
 * runs as an OpenGL 4.3 compute shader (entry points loaded at runtime via
 * wglGetProcAddress -- no SDK, no libraries). Falls back cleanly: if the
 * driver has no GL 4.3, gpu3_init returns false and the CPU solver is used. */
#ifndef GPU3D_H
#define GPU3D_H

#include "sim3d.h"
#include <stdbool.h>

/* Create shaders/buffers and upload the sim state. Call with the GL context
 * current (after gl3_init). */
bool gpu3_init(Sim3 *s);
void gpu3_shutdown(void);

/* Advance nsteps on the GPU (updates s->step). */
void gpu3_steps(Sim3 *s, int nsteps);

/* Copy the macroscopic fields back into s->rho/ux/uy/uz and fold the
 * accumulated surface force into s->fx/fy/cl/cd. steps = steps since the
 * previous readback (for force averaging). */
void gpu3_readback(Sim3 *s, int steps);

/* Mirror a geometry change: new solid mask is in s->solid; resetmask marks
 * cells that switched solid->fluid and need rest-equilibrium populations. */
void gpu3_update_geometry(Sim3 *s, const uint8_t *resetmask);

/* Reset every cell to rest equilibrium (matches sim3_reset_flow). */
void gpu3_reset_all(Sim3 *s);

#endif
