#include "gl3d.h"
#include <GL/gl.h>
#include <math.h>
#include <string.h>

static HDC   g_hdc;
static HGLRC g_ctx;
static GLuint g_tex;
static GLuint g_fontbase;
static int   g_vw = 1, g_vh = 1;

bool gl3_init(HWND hwnd)
{
    g_hdc = GetDC(hwnd);

    PIXELFORMATDESCRIPTOR pfd;
    memset(&pfd, 0, sizeof pfd);
    pfd.nSize      = sizeof pfd;
    pfd.nVersion   = 1;
    pfd.dwFlags    = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 32;
    pfd.cDepthBits = 24;
    int pf = ChoosePixelFormat(g_hdc, &pfd);
    if (!pf || !SetPixelFormat(g_hdc, pf, &pfd)) return false;
    g_ctx = wglCreateContext(g_hdc);
    if (!g_ctx || !wglMakeCurrent(g_hdc, g_ctx)) return false;

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);

    /* key light + cool fill light, world-fixed */
    const GLfloat amb[] = {0.26f, 0.26f, 0.28f, 1.f};
    const GLfloat dif[] = {0.85f, 0.85f, 0.82f, 1.f};
    glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
    const GLfloat amb1[] = {0.f, 0.f, 0.f, 1.f};
    const GLfloat dif1[] = {0.16f, 0.18f, 0.24f, 1.f};
    glLightfv(GL_LIGHT1, GL_AMBIENT, amb1);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, dif1);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    /* glossy sheen on all lit surfaces */
    const GLfloat spec[] = {0.45f, 0.45f, 0.48f, 1.f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 42.f);

    glGenTextures(1, &g_tex);
    glBindTexture(GL_TEXTURE_2D, g_tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

    /* bitmap font from the stock GUI font (ASCII 32..127) */
    g_fontbase = glGenLists(96);
    HFONT font = CreateFontA(-14, 0, 0, 0, FW_SEMIBOLD, FALSE, FALSE, FALSE,
                             ANSI_CHARSET, OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS,
                             CLEARTYPE_QUALITY, FF_DONTCARE | DEFAULT_PITCH,
                             "Segoe UI");
    HFONT old = (HFONT)SelectObject(g_hdc, font);
    wglUseFontBitmapsA(g_hdc, 32, 96, g_fontbase);
    SelectObject(g_hdc, old);

    return true;
}

void gl3_shutdown(void)
{
    if (g_ctx) { wglMakeCurrent(NULL, NULL); wglDeleteContext(g_ctx); g_ctx = NULL; }
}

void gl3_resize(int w, int h)
{
    g_vw = w > 0 ? w : 1;
    g_vh = h > 0 ? h : 1;
}

void gl3_begin_frame(const Camera3 *cam)
{
    glViewport(0, 0, g_vw, g_vh);
    glClearColor(0.02f, 0.02f, 0.045f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float aspect = (float)g_vw / (float)g_vh;
    float znear = 4.f, zfar = 4000.f;
    float top = znear * tanf(0.5f * 45.f * 3.14159265f / 180.f);
    glFrustum(-top * aspect, top * aspect, -top, top, znear, zfar);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.f, 0.f, -cam->dist);
    glRotatef(cam->pitch, 1.f, 0.f, 0.f);
    glRotatef(cam->yaw, 0.f, 1.f, 0.f);
    glTranslatef(-cam->tx, -cam->ty, -cam->tz);

    /* light directions fixed in world space */
    const GLfloat lpos[]  = {-0.45f, 0.80f, 0.40f, 0.f};
    const GLfloat lpos1[] = {0.55f, -0.45f, -0.55f, 0.f}; /* cool fill from below */
    glLightfv(GL_LIGHT0, GL_POSITION, lpos);
    glLightfv(GL_LIGHT1, GL_POSITION, lpos1);
}

void gl3_end_frame(void)
{
    SwapBuffers(g_hdc);
}

void gl3_draw_mesh(const float *nv, int nverts, float r, float g, float b)
{
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glColor3f(r, g, b);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glNormalPointer(GL_FLOAT, 6 * sizeof(float), nv);
    glVertexPointer(3, GL_FLOAT, 6 * sizeof(float), nv + 3);
    glDrawArrays(GL_TRIANGLES, 0, nverts);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_LIGHTING);
}

void gl3_draw_mesh_colored(const float *nv, const float *col3, int nverts)
{
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glNormalPointer(GL_FLOAT, 6 * sizeof(float), nv);
    glVertexPointer(3, GL_FLOAT, 6 * sizeof(float), nv + 3);
    glColorPointer(3, GL_FLOAT, 0, col3);
    glDrawArrays(GL_TRIANGLES, 0, nverts);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisable(GL_LIGHTING);
}

void gl3_draw_lines(const float *pos3, const float *col4, int nsegs, float width)
{
    if (nsegs <= 0) return;
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(width);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, pos3);
    glColorPointer(4, GL_FLOAT, 0, col4);
    glDrawArrays(GL_LINES, 0, nsegs * 2);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisable(GL_BLEND);
    glLineWidth(1.f);
}

void gl3_draw_strip(const float *pos3, const float *col4, int nverts, float width)
{
    if (nverts < 2) return;
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(width);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, pos3);
    glColorPointer(4, GL_FLOAT, 0, col4);
    glDrawArrays(GL_LINE_STRIP, 0, nverts);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisable(GL_BLEND);
    glLineWidth(1.f);
}

void gl3_draw_points(const float *pos3, const float *col4, int n, float size)
{
    if (n <= 0) return;
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE); /* additive glow */
    glDepthMask(GL_FALSE);
    glPointSize(size);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, pos3);
    glColorPointer(4, GL_FLOAT, 0, col4);
    glDrawArrays(GL_POINTS, 0, n);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
}

void gl3_draw_box(float x0, float y0, float z0, float x1, float y1, float z1)
{
    glColor4f(0.45f, 0.45f, 0.55f, 1.f);
    glBegin(GL_LINES);
    /* bottom rectangle, top rectangle, verticals */
    glVertex3f(x0,y0,z0); glVertex3f(x1,y0,z0);
    glVertex3f(x1,y0,z0); glVertex3f(x1,y0,z1);
    glVertex3f(x1,y0,z1); glVertex3f(x0,y0,z1);
    glVertex3f(x0,y0,z1); glVertex3f(x0,y0,z0);
    glVertex3f(x0,y1,z0); glVertex3f(x1,y1,z0);
    glVertex3f(x1,y1,z0); glVertex3f(x1,y1,z1);
    glVertex3f(x1,y1,z1); glVertex3f(x0,y1,z1);
    glVertex3f(x0,y1,z1); glVertex3f(x0,y1,z0);
    glVertex3f(x0,y0,z0); glVertex3f(x0,y1,z0);
    glVertex3f(x1,y0,z0); glVertex3f(x1,y1,z0);
    glVertex3f(x1,y0,z1); glVertex3f(x1,y1,z1);
    glVertex3f(x0,y0,z1); glVertex3f(x0,y1,z1);
    glEnd();
}

void gl3_draw_slice(const uint8_t *rgba, int wtex, int htex,
                    float umax, float vmax,
                    const float c0[3], const float c1[3],
                    const float c2[3], const float c3[3])
{
    glBindTexture(GL_TEXTURE_2D, g_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, wtex, htex, 0,
                 GL_RGBA, GL_UNSIGNED_BYTE, rgba);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE); /* translucent: vortex splats behind it stay visible */
    glColor4f(1.f, 1.f, 1.f, 0.88f);
    glBegin(GL_QUADS);
    glTexCoord2f(0.f,  0.f);  glVertex3fv(c0);
    glTexCoord2f(umax, 0.f);  glVertex3fv(c1);
    glTexCoord2f(umax, vmax); glVertex3fv(c2);
    glTexCoord2f(0.f,  vmax); glVertex3fv(c3);
    glEnd();
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    glDisable(GL_TEXTURE_2D);
}

void gl3_overlay_begin(int w, int h)
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, w, h, 0, -1, 1); /* top-left origin */
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glDisable(GL_DEPTH_TEST);
}

void gl3_overlay_text(int x, int y, const char *s)
{
    /* drop shadow then white text */
    glColor3f(0.f, 0.f, 0.f);
    glRasterPos2i(x + 1, y + 1);
    glListBase(g_fontbase - 32);
    glCallLists((GLsizei)strlen(s), GL_UNSIGNED_BYTE, s);
    glColor3f(1.f, 1.f, 1.f);
    glRasterPos2i(x, y);
    glCallLists((GLsizei)strlen(s), GL_UNSIGNED_BYTE, s);
}

void gl3_overlay_end(void)
{
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}
