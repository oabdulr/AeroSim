/* Software renderer: field colormaps, velocity arrows, tracer particles,
 * airfoil overlay. Draws into a 32-bit bottom-up pixel buffer that matches
 * the simulation grid one-to-one (pixel index = y * NX + x, y-up). */
#ifndef RENDER_H
#define RENDER_H

#include "sim.h"
#include "airfoil.h"

enum { VIS_PRESSURE = 0, VIS_SPEED = 1, VIS_VORTICITY = 2, VIS_COUNT = 3 };

typedef struct { float x, y; } Particle;
#define N_PART 6000

void render_field(const Sim *s, uint32_t *pix, int mode);
void render_airfoil(const Sim *s, const Airfoil *af, uint32_t *pix);
void render_arrows(const Sim *s, uint32_t *pix, int spacing);
void render_particles(const Sim *s, const Particle *p, int n, uint32_t *pix);
void render_legend(uint32_t *pix, int mode);

void particles_init(Particle *p, int n, unsigned seed);
void particles_update(const Sim *s, Particle *p, int n, float dt);

const char *vis_name(int mode);
const char *vis_range(int mode);

#endif
