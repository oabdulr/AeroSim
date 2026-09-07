#include "gpu3d.h"

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <GL/gl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* GL 4.3 subset loaded at runtime */
typedef char GLchar;
typedef ptrdiff_t GLsizeiptr;
typedef ptrdiff_t GLintptr;

#define GL_COMPUTE_SHADER              0x91B9
#define GL_SHADER_STORAGE_BUFFER       0x90D2
#define GL_SHADER_STORAGE_BARRIER_BIT  0x00002000
#define GL_BUFFER_UPDATE_BARRIER_BIT   0x00000200
#define GL_COMPILE_STATUS              0x8B81
#define GL_LINK_STATUS                 0x8B82
#define GL_STATIC_DRAW_43              0x88E4
#define GL_DYNAMIC_COPY_43             0x88EA

typedef GLuint (APIENTRY *PFN_glCreateShader)(GLenum);
typedef void   (APIENTRY *PFN_glShaderSource)(GLuint, GLsizei, const GLchar *const *, const GLint *);
typedef void   (APIENTRY *PFN_glCompileShader)(GLuint);
typedef void   (APIENTRY *PFN_glGetShaderiv)(GLuint, GLenum, GLint *);
typedef void   (APIENTRY *PFN_glGetShaderInfoLog)(GLuint, GLsizei, GLsizei *, GLchar *);
typedef GLuint (APIENTRY *PFN_glCreateProgram)(void);
typedef void   (APIENTRY *PFN_glAttachShader)(GLuint, GLuint);
typedef void   (APIENTRY *PFN_glLinkProgram)(GLuint);
typedef void   (APIENTRY *PFN_glGetProgramiv)(GLuint, GLenum, GLint *);
typedef void   (APIENTRY *PFN_glGetProgramInfoLog)(GLuint, GLsizei, GLsizei *, GLchar *);
typedef void   (APIENTRY *PFN_glUseProgram)(GLuint);
typedef void   (APIENTRY *PFN_glDeleteShader)(GLuint);
typedef void   (APIENTRY *PFN_glGenBuffers)(GLsizei, GLuint *);
typedef void   (APIENTRY *PFN_glBindBuffer)(GLenum, GLuint);
typedef void   (APIENTRY *PFN_glBufferData)(GLenum, GLsizeiptr, const void *, GLenum);
typedef void   (APIENTRY *PFN_glBufferSubData)(GLenum, GLintptr, GLsizeiptr, const void *);
typedef void   (APIENTRY *PFN_glGetBufferSubData)(GLenum, GLintptr, GLsizeiptr, void *);
typedef void   (APIENTRY *PFN_glBindBufferBase)(GLenum, GLuint, GLuint);
typedef void   (APIENTRY *PFN_glDeleteBuffers)(GLsizei, const GLuint *);
typedef void   (APIENTRY *PFN_glDispatchCompute)(GLuint, GLuint, GLuint);
typedef void   (APIENTRY *PFN_glMemoryBarrier)(GLbitfield);
typedef GLint  (APIENTRY *PFN_glGetUniformLocation)(GLuint, const GLchar *);
typedef void   (APIENTRY *PFN_glUniform1f)(GLint, GLfloat);
typedef void   (APIENTRY *PFN_glUniform1i)(GLint, GLint);

static PFN_glCreateShader       pglCreateShader;
static PFN_glShaderSource       pglShaderSource;
static PFN_glCompileShader      pglCompileShader;
static PFN_glGetShaderiv        pglGetShaderiv;
static PFN_glGetShaderInfoLog   pglGetShaderInfoLog;
static PFN_glCreateProgram      pglCreateProgram;
static PFN_glAttachShader       pglAttachShader;
static PFN_glLinkProgram        pglLinkProgram;
static PFN_glGetProgramiv       pglGetProgramiv;
static PFN_glGetProgramInfoLog  pglGetProgramInfoLog;
static PFN_glUseProgram         pglUseProgram;
static PFN_glDeleteShader       pglDeleteShader;
static PFN_glGenBuffers         pglGenBuffers;
static PFN_glBindBuffer         pglBindBuffer;
static PFN_glBufferData         pglBufferData;
static PFN_glBufferSubData      pglBufferSubData;
static PFN_glGetBufferSubData   pglGetBufferSubData;
static PFN_glBindBufferBase     pglBindBufferBase;
static PFN_glDeleteBuffers      pglDeleteBuffers;
static PFN_glDispatchCompute    pglDispatchCompute;
static PFN_glMemoryBarrier      pglMemoryBarrier;
static PFN_glGetUniformLocation pglGetUniformLocation;
static PFN_glUniform1f          pglUniform1f;
static PFN_glUniform1i          pglUniform1i;

static bool load_gl(void)
{
#define L(n) (p##n = (PFN_##n)(void *)wglGetProcAddress(#n))
    L(glCreateShader); L(glShaderSource); L(glCompileShader);
    L(glGetShaderiv); L(glGetShaderInfoLog);
    L(glCreateProgram); L(glAttachShader); L(glLinkProgram);
    L(glGetProgramiv); L(glGetProgramInfoLog); L(glUseProgram); L(glDeleteShader);
    L(glGenBuffers); L(glBindBuffer); L(glBufferData); L(glBufferSubData);
    L(glGetBufferSubData); L(glBindBufferBase); L(glDeleteBuffers);
    L(glDispatchCompute); L(glMemoryBarrier);
    L(glGetUniformLocation); L(glUniform1f); L(glUniform1i);
#undef L
    return pglCreateShader && pglDispatchCompute && pglBindBufferBase &&
           pglGetBufferSubData && pglMemoryBarrier;
}

/* shaders */

/* Shared GLSL prologue: constants + equilibrium. */
static const char *SRC_COMMON =
"const ivec3 E[19]=ivec3[19](ivec3(0,0,0),\n"
"  ivec3(1,0,0),ivec3(-1,0,0),ivec3(0,1,0),ivec3(0,-1,0),ivec3(0,0,1),ivec3(0,0,-1),\n"
"  ivec3(1,1,0),ivec3(-1,-1,0),ivec3(1,-1,0),ivec3(-1,1,0),\n"
"  ivec3(1,0,1),ivec3(-1,0,-1),ivec3(1,0,-1),ivec3(-1,0,1),\n"
"  ivec3(0,1,1),ivec3(0,-1,-1),ivec3(0,1,-1),ivec3(0,-1,1));\n"
"const float W[19]=float[19](1.0/3.0,\n"
"  1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,\n"
"  1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,\n"
"  1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0);\n"
"const int OPP[19]=int[19](0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17);\n"
"float feq(int q,float rho,vec3 u){\n"
"  float eu=dot(vec3(E[q]),u); float u2=dot(u,u);\n"
"  return W[q]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*u2);\n"
"}\n";

/* Fused stream(pull)-collide with LES, absorbing layers, bounce-back and
 * fixed-point force accumulation. One thread per cell. */
static const char *SRC_STEP =
"layout(local_size_x=128) in;\n"
"layout(std430,binding=0) readonly  buffer FA { float fa[]; };\n"
"layout(std430,binding=1) writeonly buffer FB { float fb[]; };\n"
"layout(std430,binding=2) readonly  buffer FL { uint  flg[]; };\n"
"layout(std430,binding=3) readonly  buffer SG { float sig[]; };\n"
"layout(std430,binding=4) buffer MC { vec4 mac[]; };\n"
"layout(std430,binding=5) buffer FR { int  frc[]; };\n"
"uniform float uin; uniform float tau0; uniform float smag2;\n"
"void main(){\n"
"  uint gid=gl_GlobalInvocationID.x; if(gid>=uint(NC)) return;\n"
"  int i=int(gid);\n"
"  if((flg[i]&1u)!=0u) return;\n"
"  int x=i%NX; int y=(i/NX)%NY; int z=i/(NX*NY);\n"
"  float f[19];\n"
"  for(int q=0;q<19;q++){\n"
"    int xs=x-E[q].x;\n"
"    if(xs<0||xs>=NX){ f[q]=fa[q*NC+i]; continue; }\n"
"    int ys=y-E[q].y; if(ys<0) ys+=NY; if(ys>=NY) ys-=NY;\n"
"    int zs=z-E[q].z; if(zs<0) zs+=NZ; if(zs>=NZ) zs-=NZ;\n"
"    int j=(zs*NY+ys)*NX+xs;\n"
"    if((flg[j]&1u)!=0u){\n"
"      float fv=fa[OPP[q]*NC+i];\n"
"      f[q]=fv;\n"
"      atomicAdd(frc[0],int(-2.0*float(E[q].x)*fv*1.0e6));\n"
"      atomicAdd(frc[1],int(-2.0*float(E[q].y)*fv*1.0e6));\n"
"      atomicAdd(frc[2],int(-2.0*float(E[q].z)*fv*1.0e6));\n"
"    } else f[q]=fa[q*NC+j];\n"
"  }\n"
"  float rho=0.0; vec3 m=vec3(0.0);\n"
"  for(int q=0;q<19;q++){ rho+=f[q]; m+=vec3(E[q])*f[q]; }\n"
"  vec3 u=m/rho;\n"
"  float um2=dot(u,u); if(um2>0.15) u*=sqrt(0.15/um2);\n"
"  float fe[19];\n"
"  float pxx=0.0,pyy=0.0,pzz=0.0,pxy=0.0,pxz=0.0,pyz=0.0;\n"
"  for(int q=0;q<19;q++){\n"
"    fe[q]=feq(q,rho,u);\n"
"    float fn=f[q]-fe[q]; vec3 e=vec3(E[q]);\n"
"    pxx+=e.x*e.x*fn; pyy+=e.y*e.y*fn; pzz+=e.z*e.z*fn;\n"
"    pxy+=e.x*e.y*fn; pxz+=e.x*e.z*fn; pyz+=e.y*e.z*fn;\n"
"  }\n"
"  float qn=sqrt(pxx*pxx+pyy*pyy+pzz*pzz+2.0*(pxy*pxy+pxz*pxz+pyz*pyz));\n"
"  float taue=0.5*(tau0+sqrt(tau0*tau0+25.4558441*smag2*qn/rho));\n"
"  float om=1.0/taue;\n"
"  float sg=sig[i];\n"
"  for(int q=0;q<19;q++){\n"
"    float v=f[q]-om*(f[q]-fe[q]);\n"
"    if(sg>0.0) v-=sg*(v-feq(q,1.0,vec3(uin,0.0,0.0)));\n"
"    fb[q*NC+i]=v;\n"
"  }\n"
"  mac[i]=vec4(rho,u);\n"
"}\n";

/* Inlet/outlet planes; one thread per (y,z) column pair. */
static const char *SRC_BC =
"layout(local_size_x=128) in;\n"
"layout(std430,binding=1) buffer FB { float fb[]; };\n"
"layout(std430,binding=4) buffer MC { vec4 mac[]; };\n"
"uniform float uin;\n"
"void main(){\n"
"  uint gid=gl_GlobalInvocationID.x; if(gid>=uint(NY*NZ)) return;\n"
"  int y=int(gid)%NY; int z=int(gid)/NY;\n"
"  int i0=(z*NY+y)*NX;\n"
"  float r=clamp(mac[i0+1].x,0.5,2.0);\n"
"  for(int q=0;q<19;q++) fb[q*NC+i0]=feq(q,r,vec3(uin,0.0,0.0));\n"
"  mac[i0]=vec4(r,uin,0.0,0.0);\n"
"  int i1=i0+NX-1;\n"
"  for(int q=0;q<19;q++) fb[q*NC+i1]=fb[q*NC+i1-1];\n"
"  mac[i1]=mac[i1-1];\n"
"}\n";

/* Rest-equilibrium reset: mode 0 = only cells flagged with bit 1,
 * mode 1 = every cell. */
static const char *SRC_RESET =
"layout(local_size_x=128) in;\n"
"layout(std430,binding=0) buffer FA { float fa[]; };\n"
"layout(std430,binding=1) buffer FB { float fb[]; };\n"
"layout(std430,binding=2) readonly buffer FL { uint flg[]; };\n"
"layout(std430,binding=4) buffer MC { vec4 mac[]; };\n"
"uniform int mode;\n"
"void main(){\n"
"  uint gid=gl_GlobalInvocationID.x; if(gid>=uint(NC)) return;\n"
"  int i=int(gid);\n"
"  if(mode==0 && (flg[i]&2u)==0u) return;\n"
"  for(int q=0;q<19;q++){ fa[q*NC+i]=W[q]; fb[q*NC+i]=W[q]; }\n"
"  mac[i]=vec4(1.0,0.0,0.0,0.0);\n"
"}\n";

/* backend state */
static GLuint g_progStep, g_progBC, g_progReset;
static GLint  g_uStepUin, g_uStepTau, g_uStepSmag, g_uBCUin, g_uResetMode;
static GLuint g_bufF[2], g_bufFlags, g_bufSigma, g_bufMacro, g_bufForce;
static int    g_src;                   /* which of bufF is the source */
static float *g_stage;                 /* readback staging, NCELLS3 vec4 */
static unsigned *g_flagstage;
static long long g_forceSteps;
static bool g_ready;

static GLuint compile_program(const char *body, const char *name)
{
    char header[256];
    snprintf(header, sizeof header,
             "#version 430\n#define NX %d\n#define NY %d\n#define NZ %d\n#define NC %d\n",
             NX3, NY3, NZ3, NCELLS3);
    const char *parts[3] = { header, SRC_COMMON, body };

    GLuint sh = pglCreateShader(GL_COMPUTE_SHADER);
    pglShaderSource(sh, 3, parts, NULL);
    pglCompileShader(sh);
    GLint ok = 0;
    pglGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[2048]; GLsizei n = 0;
        pglGetShaderInfoLog(sh, sizeof log, &n, log);
        char msg[2300];
        snprintf(msg, sizeof msg, "GPU shader '%s' failed to compile:\n%.2000s", name, log);
        MessageBoxA(NULL, msg, "AeroSim3D", MB_ICONWARNING);
        return 0;
    }
    GLuint pr = pglCreateProgram();
    pglAttachShader(pr, sh);
    pglLinkProgram(pr);
    pglDeleteShader(sh);
    pglGetProgramiv(pr, GL_LINK_STATUS, &ok);
    if (!ok) return 0;
    return pr;
}

static void upload_flags(const Sim3 *s, const uint8_t *resetmask)
{
    for (int i = 0; i < NCELLS3; ++i)
        g_flagstage[i] = (unsigned)(s->solid[i] ? 1u : 0u)
                       | (resetmask && resetmask[i] ? 2u : 0u);
    pglBindBuffer(GL_SHADER_STORAGE_BUFFER, g_bufFlags);
    pglBufferSubData(GL_SHADER_STORAGE_BUFFER, 0,
                     sizeof(unsigned) * NCELLS3, g_flagstage);
}

bool gpu3_init(Sim3 *s)
{
    const char *ver = (const char *)glGetString(GL_VERSION);
    if (!ver || !load_gl()) return false;
    int maj = 0, min = 0;
    sscanf(ver, "%d.%d", &maj, &min);
    if (maj < 4 || (maj == 4 && min < 3)) return false;

    g_progStep  = compile_program(SRC_STEP,  "step");
    g_progBC    = compile_program(SRC_BC,    "boundary");
    g_progReset = compile_program(SRC_RESET, "reset");
    if (!g_progStep || !g_progBC || !g_progReset) return false;
    g_uStepUin  = pglGetUniformLocation(g_progStep, "uin");
    g_uStepTau  = pglGetUniformLocation(g_progStep, "tau0");
    g_uStepSmag = pglGetUniformLocation(g_progStep, "smag2");
    g_uBCUin    = pglGetUniformLocation(g_progBC, "uin");
    g_uResetMode = pglGetUniformLocation(g_progReset, "mode");

    g_stage     = (float *)malloc(sizeof(float) * 4 * NCELLS3);
    g_flagstage = (unsigned *)malloc(sizeof(unsigned) * NCELLS3);
    if (!g_stage || !g_flagstage) return false;

    pglGenBuffers(2, g_bufF);
    pglGenBuffers(1, &g_bufFlags);
    pglGenBuffers(1, &g_bufSigma);
    pglGenBuffers(1, &g_bufMacro);
    pglGenBuffers(1, &g_bufForce);

    for (int k = 0; k < 2; ++k) {
        pglBindBuffer(GL_SHADER_STORAGE_BUFFER, g_bufF[k]);
        pglBufferData(GL_SHADER_STORAGE_BUFFER,
                      (GLsizeiptr)sizeof(float) * Q3 * NCELLS3,
                      k == 0 ? s->f : s->fnew, GL_DYNAMIC_COPY_43);
    }
    pglBindBuffer(GL_SHADER_STORAGE_BUFFER, g_bufFlags);
    pglBufferData(GL_SHADER_STORAGE_BUFFER, (GLsizeiptr)sizeof(unsigned) * NCELLS3,
                  NULL, GL_DYNAMIC_COPY_43);
    upload_flags(s, NULL);

    pglBindBuffer(GL_SHADER_STORAGE_BUFFER, g_bufSigma);
    pglBufferData(GL_SHADER_STORAGE_BUFFER, (GLsizeiptr)sizeof(float) * NCELLS3,
                  s->sigma, GL_STATIC_DRAW_43);

    for (int i = 0; i < NCELLS3; ++i) {
        g_stage[i * 4 + 0] = s->rho[i];
        g_stage[i * 4 + 1] = s->ux[i];
        g_stage[i * 4 + 2] = s->uy[i];
        g_stage[i * 4 + 3] = s->uz[i];
    }
    pglBindBuffer(GL_SHADER_STORAGE_BUFFER, g_bufMacro);
    pglBufferData(GL_SHADER_STORAGE_BUFFER, (GLsizeiptr)sizeof(float) * 4 * NCELLS3,
                  g_stage, GL_DYNAMIC_COPY_43);

    const int zero[4] = {0, 0, 0, 0};
    pglBindBuffer(GL_SHADER_STORAGE_BUFFER, g_bufForce);
    pglBufferData(GL_SHADER_STORAGE_BUFFER, sizeof zero, zero, GL_DYNAMIC_COPY_43);

    g_src = 0;
    g_forceSteps = 0;
    g_ready = true;
    return true;
}

void gpu3_shutdown(void)
{
    if (!g_ready) return;
    pglDeleteBuffers(2, g_bufF);
    pglDeleteBuffers(1, &g_bufFlags);
    pglDeleteBuffers(1, &g_bufSigma);
    pglDeleteBuffers(1, &g_bufMacro);
    pglDeleteBuffers(1, &g_bufForce);
    free(g_stage); free(g_flagstage);
    g_ready = false;
}

void gpu3_steps(Sim3 *s, int nsteps)
{
    const GLuint groupsCell = (GLuint)((NCELLS3 + 127) / 128);
    const GLuint groupsBC   = (GLuint)((NY3 * NZ3 + 127) / 128);

    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, g_bufFlags);
    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, g_bufSigma);
    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, g_bufMacro);
    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, g_bufForce);

    for (int k = 0; k < nsteps; ++k) {
        float uin = sim3_inlet_speed(s);

        pglUseProgram(g_progStep);
        pglUniform1f(g_uStepUin, uin);
        pglUniform1f(g_uStepTau, s->tau0);
        pglUniform1f(g_uStepSmag, s->smag2);
        pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, g_bufF[g_src]);
        pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, g_bufF[1 - g_src]);
        pglDispatchCompute(groupsCell, 1, 1);
        pglMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        pglUseProgram(g_progBC);
        pglUniform1f(g_uBCUin, uin);
        pglDispatchCompute(groupsBC, 1, 1);
        pglMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        g_src = 1 - g_src;
        s->step++;
        g_forceSteps++;
    }
}

void gpu3_readback(Sim3 *s, int steps)
{
    (void)steps;
    pglBindBuffer(GL_SHADER_STORAGE_BUFFER, g_bufMacro);
    pglGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0,
                        (GLsizeiptr)sizeof(float) * 4 * NCELLS3, g_stage);
    for (int i = 0; i < NCELLS3; ++i) {
        s->rho[i] = g_stage[i * 4 + 0];
        s->ux[i]  = g_stage[i * 4 + 1];
        s->uy[i]  = g_stage[i * 4 + 2];
        s->uz[i]  = g_stage[i * 4 + 3];
    }

    if (g_forceSteps > 0) {
        int acc[3] = {0, 0, 0};
        pglBindBuffer(GL_SHADER_STORAGE_BUFFER, g_bufForce);
        pglGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof acc, acc);
        const int zero[4] = {0, 0, 0, 0};
        pglBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof zero, zero);

        double inv = 1.0 / (1.0e6 * (double)g_forceSteps);
        s->fx = acc[0] * inv;
        s->fy = acc[1] * inv;
        s->fz = acc[2] * inv;
        g_forceSteps = 0;

        float uin  = sim3_inlet_speed(s);
        float uref = uin > 0.02f ? uin : 0.02f;
        double qd  = 0.5 * (double)uref * uref * s->refS;
        double a   = 0.02; /* per-frame smoothing */
        s->cl += a * (s->fy / qd - s->cl);
        s->cd += a * (s->fx / qd - s->cd);
    }
}

void gpu3_update_geometry(Sim3 *s, const uint8_t *resetmask)
{
    if (!g_ready) return;
    upload_flags(s, resetmask);
    pglUseProgram(g_progReset);
    pglUniform1i(g_uResetMode, 0);
    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, g_bufF[g_src]);
    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, g_bufF[1 - g_src]);
    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, g_bufFlags);
    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, g_bufMacro);
    pglDispatchCompute((GLuint)((NCELLS3 + 127) / 128), 1, 1);
    pglMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_BUFFER_UPDATE_BARRIER_BIT);
    upload_flags(s, NULL); /* clear the reset bits */
}

void gpu3_reset_all(Sim3 *s)
{
    (void)s;
    if (!g_ready) return;
    pglUseProgram(g_progReset);
    pglUniform1i(g_uResetMode, 1);
    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, g_bufF[g_src]);
    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, g_bufF[1 - g_src]);
    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, g_bufFlags);
    pglBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, g_bufMacro);
    pglDispatchCompute((GLuint)((NCELLS3 + 127) / 128), 1, 1);
    pglMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_BUFFER_UPDATE_BARRIER_BIT);

    const int zero[4] = {0, 0, 0, 0};
    pglBindBuffer(GL_SHADER_STORAGE_BUFFER, g_bufForce);
    pglBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof zero, zero);
    g_forceSteps = 0;
}
