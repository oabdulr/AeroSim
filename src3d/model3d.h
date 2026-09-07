/* Triangle-mesh model import (binary/ASCII STL and Wavefront OBJ) with
 * voxelization into the simulation's solid mask. Lets any watertight 3D
 * model be placed in the flow instead of the built-in wing. */
#ifndef MODEL3D_H
#define MODEL3D_H

#include "sim3d.h"
#include <stdbool.h>

/* Load a model file. On failure returns false and writes a reason to err. */
bool model3_load(const char *path, char *err, int errlen);

bool        model3_active(void);
void        model3_clear(void);       /* back to the built-in wing        */
void        model3_cycle_axes(void);  /* fix models authored z-up etc.    */
const char *model3_name(void);
double      model3_refS(void);        /* frontal-area proxy for CL/CD     */
int         model3_tris(void);

/* Voxelize the model (yawed about the vertical axis, then pitched nose-up)
 * into s->solid, refreshing changed cells and the careful mask. */
void model3_rasterize(Sim3 *s, float pitch_deg, float yaw_deg);

/* Render mesh, interleaved [n xyz | v xyz], 3 verts/tri. Returns verts. */
int model3_mesh(float *out, int maxverts, float pitch_deg, float yaw_deg);

#endif
