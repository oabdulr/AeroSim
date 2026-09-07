/* AeroSim -- interactive 2D LBM wind tunnel around a NACA airfoil.
 * Pure Win32 + GDI, no external dependencies. See README.md for physics. */
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "sim.h"
#include "airfoil.h"
#include "render.h"

static Sim      g_sim;
static Airfoil  g_af;
static Particle g_part[N_PART];
static uint32_t *g_pix;
static HWND     g_hwnd;
static int      g_cw = 1280, g_ch = 640;   /* client size */

static int   g_mode      = VIS_PRESSURE;
static BOOL  g_arrows    = TRUE;
static BOOL  g_particles = TRUE;
static BOOL  g_paused    = FALSE;
static BOOL  g_help      = TRUE;
static BOOL  g_run       = TRUE;
static int   g_spf       = 12;      /* sim steps per rendered frame */
static float g_flap_mem  = 20.f;    /* remembered deflection for F toggle */
static double g_fps, g_mlups;

static void rebuild_airfoil(void)
{
    airfoil_build(&g_af);
    airfoil_rasterize(&g_af, &g_sim);
}

static LRESULT CALLBACK wndproc(HWND h, UINT m, WPARAM w, LPARAM l)
{
    switch (m) {
    case WM_SIZE:
        g_cw = LOWORD(l);
        g_ch = HIWORD(l);
        return 0;
    case WM_ERASEBKGND:
        return 1;
    case WM_PAINT: {
        PAINTSTRUCT ps;
        BeginPaint(h, &ps);
        EndPaint(h, &ps);
        return 0;
    }
    case WM_CLOSE:
    case WM_DESTROY:
        g_run = FALSE;
        PostQuitMessage(0);
        return 0;
    case WM_KEYDOWN:
        switch (w) {
        case VK_ESCAPE: g_run = FALSE; PostQuitMessage(0); break;
        case VK_UP:
            g_af.aoa_deg += 0.5f;
            if (g_af.aoa_deg > 25.f) g_af.aoa_deg = 25.f;
            rebuild_airfoil();
            break;
        case VK_DOWN:
            g_af.aoa_deg -= 0.5f;
            if (g_af.aoa_deg < -25.f) g_af.aoa_deg = -25.f;
            rebuild_airfoil();
            break;
        case 'F':
            if (g_af.flap_deg != 0.f) {
                g_flap_mem = g_af.flap_deg;
                g_af.flap_deg = 0.f;
            } else {
                g_af.flap_deg = g_flap_mem;
            }
            rebuild_airfoil();
            break;
        case VK_OEM_4: /* [ */
            g_af.flap_deg -= 5.f;
            if (g_af.flap_deg < -20.f) g_af.flap_deg = -20.f;
            rebuild_airfoil();
            break;
        case VK_OEM_6: /* ] */
            g_af.flap_deg += 5.f;
            if (g_af.flap_deg > 40.f) g_af.flap_deg = 40.f;
            rebuild_airfoil();
            break;
        case '1': g_mode = VIS_PRESSURE;  break;
        case '2': g_mode = VIS_SPEED;     break;
        case '3': g_mode = VIS_VORTICITY; break;
        case 'A': g_arrows = !g_arrows;       break;
        case 'P': g_particles = !g_particles; break;
        case 'H': g_help = !g_help;           break;
        case 'R': sim_reset_flow(&g_sim);     break;
        case VK_SPACE: g_paused = !g_paused;  break;
        case VK_OEM_COMMA: { /* , = lower Reynolds */
            double re = g_sim.reynolds / 1.5;
            if (re < 500.0) re = 500.0;
            sim_set_reynolds(&g_sim, re);
            break;
        }
        case VK_OEM_PERIOD: { /* . = raise Reynolds */
            double re = g_sim.reynolds * 1.5;
            if (re > 100000.0) re = 100000.0;
            sim_set_reynolds(&g_sim, re);
            break;
        }
        case VK_OEM_MINUS: case VK_SUBTRACT:
            g_spf -= 2; if (g_spf < 1) g_spf = 1;
            break;
        case VK_OEM_PLUS: case VK_ADD:
            g_spf += 2; if (g_spf > 48) g_spf = 48;
            break;
        }
        return 0;
    }
    return DefWindowProcA(h, m, w, l);
}

static void hud_line(HDC hdc, int x, int y, const char *txt)
{
    SetTextColor(hdc, RGB(0, 0, 0));
    TextOutA(hdc, x + 1, y + 1, txt, (int)strlen(txt));
    SetTextColor(hdc, RGB(255, 255, 255));
    TextOutA(hdc, x, y, txt, (int)strlen(txt));
}

static void present(HDC hdc, int steps_this_frame)
{
    (void)steps_this_frame;
    if (g_cw < 8 || g_ch < 8) return;

    BITMAPINFO bmi;
    memset(&bmi, 0, sizeof bmi);
    bmi.bmiHeader.biSize        = sizeof(BITMAPINFOHEADER);
    bmi.bmiHeader.biWidth       = NX;
    bmi.bmiHeader.biHeight      = NY; /* positive = bottom-up = grid y-up */
    bmi.bmiHeader.biPlanes      = 1;
    bmi.bmiHeader.biBitCount    = 32;
    bmi.bmiHeader.biCompression = BI_RGB;

    /* Letterbox to preserve the 2:1 tunnel aspect. */
    double sc = (double)g_cw / NX;
    double s2 = (double)g_ch / NY;
    if (s2 < sc) sc = s2;
    int dw = (int)(NX * sc), dh = (int)(NY * sc);
    int ox = (g_cw - dw) / 2, oy = (g_ch - dh) / 2;

    SetStretchBltMode(hdc, HALFTONE);
    SetBrushOrgEx(hdc, 0, 0, NULL);
    StretchDIBits(hdc, ox, oy, dw, dh, 0, 0, NX, NY,
                  g_pix, &bmi, DIB_RGB_COLORS, SRCCOPY);

    /* Black bars around the letterboxed image (unconditional, so rounding
     * never leaves an unpainted 1px strip at the right/bottom edges). */
    HBRUSH black = (HBRUSH)GetStockObject(BLACK_BRUSH);
    RECT r;
    if (oy > 0)          { r = (RECT){0, 0, g_cw, oy};          FillRect(hdc, &r, black); }
    if (oy + dh < g_ch)  { r = (RECT){0, oy + dh, g_cw, g_ch};  FillRect(hdc, &r, black); }
    if (ox > 0)          { r = (RECT){0, 0, ox, g_ch};          FillRect(hdc, &r, black); }
    if (ox + dw < g_cw)  { r = (RECT){ox + dw, 0, g_cw, g_ch};  FillRect(hdc, &r, black); }

    /* HUD */
    SetBkMode(hdc, TRANSPARENT);
    char buf[320];
    int  tx = ox + 10, ty = oy + 8;

    double tstar = (double)g_sim.step * g_sim.u0 / g_sim.chord;
    snprintf(buf, sizeof buf,
             "NACA %02d%02d  AoA %+.1f\xB0   Flap %+.0f\xB0   Re %.0f   t* %.1f   view: %s",
             (int)(g_af.camber * 100.f + 0.5f) * 10 + (int)(g_af.camber_pos * 10.f + 0.5f),
             (int)(g_af.thickness * 100.f + 0.5f),
             (double)g_af.aoa_deg, (double)g_af.flap_deg,
             g_sim.reynolds, tstar, vis_name(g_mode));
    hud_line(hdc, tx, ty, buf); ty += 17;

    snprintf(buf, sizeof buf,
             "Cl %+.3f   Cd %+.4f   L/D %+.1f   |   %d steps/frame   %.0f fps   %.0f MLUPS%s%s",
             g_sim.cl, g_sim.cd,
             fabs(g_sim.cd) > 1e-6 ? g_sim.cl / g_sim.cd : 0.0,
             g_spf, g_fps, g_mlups,
             g_paused ? "   [PAUSED]" : "",
             g_sim.step < g_sim.ramp_steps ? "   [flow spinning up...]" : "");
    hud_line(hdc, tx, ty, buf); ty += 17;

    if (g_help) {
        hud_line(hdc, tx, ty,
                 "UP/DOWN angle of attack   F flap toggle   [ ] flap angle   1 pressure  2 velocity  3 vorticity");
        ty += 17;
        hud_line(hdc, tx, ty,
                 "A arrows   P particles   , . Reynolds   +/- sim speed   SPACE pause   R reset flow   H hide help");
        ty += 17;
    }

    /* Legend label, anchored just above the color bar (grid y0=14, h=12). */
    int ly = oy + dh - (int)((14 + 12 + 18) * sc / 1.0);
    hud_line(hdc, ox + (int)(14 * sc), ly, vis_range(g_mode));
}

int WINAPI WinMain(HINSTANCE hi, HINSTANCE hp, LPSTR cmd, int show)
{
    (void)hp; (void)cmd; (void)show;
    SetProcessDPIAware();

    WNDCLASSA wc;
    memset(&wc, 0, sizeof wc);
    wc.style         = CS_OWNDC;
    wc.lpfnWndProc   = wndproc;
    wc.hInstance     = hi;
    wc.hCursor       = LoadCursor(NULL, IDC_ARROW);
    wc.lpszClassName = "AeroSimWnd";
    RegisterClassA(&wc);

    RECT wr = {0, 0, g_cw, g_ch};
    AdjustWindowRect(&wr, WS_OVERLAPPEDWINDOW, FALSE);
    g_hwnd = CreateWindowExA(0, "AeroSimWnd",
                             "AeroSim - 2D open-air flow (D2Q9 LBM + Smagorinsky LES)",
                             WS_OVERLAPPEDWINDOW | WS_VISIBLE,
                             CW_USEDEFAULT, CW_USEDEFAULT,
                             wr.right - wr.left, wr.bottom - wr.top,
                             NULL, NULL, hi, NULL);
    if (!g_hwnd) return 1;

    if (!sim_init(&g_sim)) {
        MessageBoxA(NULL, "Out of memory allocating the simulation grid.",
                    "AeroSim", MB_ICONERROR);
        return 1;
    }
    g_pix = (uint32_t *)malloc(sizeof(uint32_t) * NCELLS);
    if (!g_pix) return 1;

    memset(&g_af, 0, sizeof g_af);
    g_af.camber     = 0.02f;   /* NACA 2412 */
    g_af.camber_pos = 0.4f;
    g_af.thickness  = 0.12f;
    g_af.aoa_deg    = 6.f;
    g_af.flap_deg   = 0.f;
    g_af.flap_hinge = 0.70f;
    g_af.chord      = g_sim.chord;
    g_af.cx         = NX * 0.35f;
    g_af.cy         = NY * 0.5f;
    rebuild_airfoil();
    particles_init(g_part, N_PART, 20260712u);

    LARGE_INTEGER freq, t0, t1;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&t0);
    HDC hdc = GetDC(g_hwnd);

    MSG msg;
    while (g_run) {
        while (PeekMessageA(&msg, NULL, 0, 0, PM_REMOVE)) {
            if (msg.message == WM_QUIT) g_run = FALSE;
            TranslateMessage(&msg);
            DispatchMessageA(&msg);
        }
        if (!g_run) break;

        int steps = 0;
        if (!g_paused) {
            for (int k = 0; k < g_spf; ++k) sim_step(&g_sim);
            steps = g_spf;
            if (g_particles) particles_update(&g_sim, g_part, N_PART, (float)steps);
        }

        /* Watchdog: if the field ever goes non-finite, restart the flow
         * instead of filling the screen with NaN garbage. The NaN/Inf test
         * inspects the exponent bits directly so it survives /fp:fast and
         * -ffast-math, which are allowed to assume floats are finite. */
        float probe = g_sim.rho[(NY / 2) * NX + NX / 2];
        uint32_t pbits;
        memcpy(&pbits, &probe, sizeof pbits);
        int nonfinite = ((pbits >> 23) & 0xFFu) == 0xFFu;
        if (nonfinite || probe < 0.1f || probe > 10.f) sim_reset_flow(&g_sim);

        render_field(&g_sim, g_pix, g_mode);
        render_airfoil(&g_sim, &g_af, g_pix);
        if (g_arrows) render_arrows(&g_sim, g_pix, 28);
        if (g_particles) render_particles(&g_sim, g_part, N_PART, g_pix);
        render_legend(g_pix, g_mode);
        present(hdc, steps);

        QueryPerformanceCounter(&t1);
        double dt = (double)(t1.QuadPart - t0.QuadPart) / (double)freq.QuadPart;
        t0 = t1;
        if (dt > 1e-6) {
            g_fps   += 0.1 * (1.0 / dt - g_fps);
            g_mlups += 0.1 * ((double)steps * NCELLS / dt / 1e6 - g_mlups);
        }
        if (g_paused) Sleep(15);
    }

    ReleaseDC(g_hwnd, hdc);
    sim_free(&g_sim);
    free(g_pix);
    return 0;
}
