/* Minimal fixed-pipeline OpenGL 1.1 renderer (opengl32.dll ships with
 * Windows, so this adds no dependency). Camera orbit, lit meshes, point
 * clouds, one dynamic 256x128 RGBA texture for the field slice, and
 * bitmap-font overlay text. */
#ifndef GL3D_H
#define GL3D_H

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <stdint.h>
#include <stdbool.h>

typedef struct {
    float yaw, pitch, dist;
    float tx, ty, tz; /* orbit target */
} Camera3;

bool gl3_init(HWND hwnd);
void gl3_shutdown(void);
void gl3_resize(int w, int h);

void gl3_begin_frame(const Camera3 *cam);
void gl3_end_frame(void);

void gl3_draw_mesh(const float *nv_interleaved, int nverts,
                   float r, float g, float b);
/* Same mesh layout plus one RGB per vertex (surface-pressure coloring). */
void gl3_draw_mesh_colored(const float *nv_interleaved, const float *col3,
                           int nverts);
void gl3_draw_points(const float *pos3, const float *col4, int n, float size);
/* n line segments: pos3 has 2 vertices per segment, col4 one RGBA per vertex. */
void gl3_draw_lines(const float *pos3, const float *col4, int nsegs, float width);
/* One polyline with per-vertex RGBA colors. */
void gl3_draw_strip(const float *pos3, const float *col4, int nverts, float width);
void gl3_draw_box(float x0, float y0, float z0, float x1, float y1, float z1);

/* Upload a wtex x htex RGBA byte image and draw it as a quad with corners
 * c0..c3 (counter-clockwise); u/v extents give the used sub-rectangle. */
void gl3_draw_slice(const uint8_t *rgba, int wtex, int htex,
                    float umax, float vmax,
                    const float c0[3], const float c1[3],
                    const float c2[3], const float c3[3]);

void gl3_overlay_begin(int w, int h);
void gl3_overlay_text(int x, int y, const char *s);
void gl3_overlay_end(void);

#endif
