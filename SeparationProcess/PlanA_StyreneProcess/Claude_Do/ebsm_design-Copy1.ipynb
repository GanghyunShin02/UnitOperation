"""
Ethylbenzene (EB) / styrene monomer (SM) distillation column design at 1 atm.

Sequence (lecture notes): VLE -> material balance -> minimum reflux -> minimum
stages (Fenske) -> McCabe-Thiele -> stage efficiency (O'Connell / Murphree)
-> packed-column alternative (NTU, HTU).

Requires: numpy, scipy, matplotlib.  Run:  python ebsm_design.py
Figures are written to ./images/.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

# ----------------------------------------------------------------------------
# 0. Design basis (edit here)
# ----------------------------------------------------------------------------
P = 760.0            # column pressure, mmHg (uniform; pressure drop neglected)
F = 100.0            # feed rate, kmol/h (calculation basis)
xF, xD, xB = 0.50, 0.99, 0.01   # mole fraction EB (light component)
q = 1.0              # saturated-liquid feed
R_FACTOR = 1.2       # R_D = R_FACTOR * R_Dm

# Antoine constants: log10(P[mmHg]) = A - B / (T[degC] + C)
ANT = {"EB": (6.95719, 1424.255, 213.206), "SM": (6.92409, 1420.0, 206.0)}
MW = {"EB": 106.165, "SM": 104.149}

# Pure-component properties at 140 degC, 1 atm (thermo library correlations)
RHO_L = {"EB": 756.5, "SM": 789.9}        # liquid density, kg/m3
MU_L = {"EB": 0.227, "SM": 0.234}         # liquid viscosity, cP
MU_V = {"EB": 8.91e-6, "SM": 8.82e-6}     # vapor viscosity, Pa.s
VB = {"EB": 139.6, "SM": 132.8}           # liquid molar volume at Tb, cm3/mol
# Fuller-Schettler-Giddings atomic diffusion volumes: C 16.5, H 1.98, aromatic ring -20.2
FSG_V = {"EB": 8 * 16.5 + 10 * 1.98 - 20.2, "SM": 8 * 16.5 + 8 * 1.98 - 20.2}

OUT = Path(__file__).parent / "images"
OUT.mkdir(exist_ok=True)

# ----------------------------------------------------------------------------
# 1. Vapor-liquid equilibrium: Antoine + Raoult's law
# ----------------------------------------------------------------------------
def psat(c, T):
    A, B, C = ANT[c]
    return 10 ** (A - B / (T + C))


def T_bubble(x):
    return brentq(lambda T: x * psat("EB", T) + (1 - x) * psat("SM", T) - P, 100, 180)


def T_dew(y):
    return brentq(lambda T: y * P / psat("EB", T) + (1 - y) * P / psat("SM", T) - 1, 100, 180)


def y_eq(x):
    """Vapor composition in equilibrium with liquid x (bubble-point calculation)."""
    return x * psat("EB", T_bubble(x)) / P


def x_eq(y):
    """Liquid composition in equilibrium with vapor y (dew-point calculation)."""
    return y * P / psat("EB", T_dew(y))


def alpha(x):
    T = T_bubble(x)
    return psat("EB", T) / psat("SM", T)


xg = np.linspace(0, 1, 2001)
yg = np.array([y_eq(x) for x in xg])
Tg = np.array([T_bubble(x) for x in xg])
ye_spline = CubicSpline(xg, yg)          # used only for the slope m = dy_e/dx
m_slope = ye_spline.derivative()

aD, aF, aB = alpha(xD), alpha(xF), alpha(xB)
a_mean = (aD * aF * aB) ** (1 / 3)

# ----------------------------------------------------------------------------
# 2. Overall material balance
# ----------------------------------------------------------------------------
D = F * (xF - xB) / (xD - xB)
B = F - D

# ----------------------------------------------------------------------------
# 3. Minimum reflux ratio (q = 1: feed line is vertical, x' = xF)
# ----------------------------------------------------------------------------
x_p = xF
y_p = y_eq(x_p)
R_min = (xD - y_p) / (y_p - x_p)

# ----------------------------------------------------------------------------
# 4. Minimum number of stages (total reflux)
# ----------------------------------------------------------------------------
N_min_fenske = np.log(xD * (1 - xB) / (xB * (1 - xD))) / np.log(a_mean) - 1  # plates, reboiler excluded


def step_total_reflux():
    pts, x, y, n = [(xD, xD)], xD, xD, 0
    while True:
        xn = x_eq(y)
        n += 1
        pts.append((xn, y))
        if xn <= xB:
            return n - 1 + (x - xB) / (x - xn), pts
        x = y = xn
        pts.append((x, y))


N_min_graph, pts_total = step_total_reflux()   # stages incl. reboiler

# ----------------------------------------------------------------------------
# 5. McCabe-Thiele at the operating reflux ratio
# ----------------------------------------------------------------------------
def lines(R):
    """Operating lines for q = 1. Returns rect(x), strip(x), their inverses, and the intersection."""
    xi = xF
    yi = R / (R + 1) * xi + xD / (R + 1)
    rect = lambda x: R / (R + 1) * x + xD / (R + 1)
    strip = lambda x: xB + (x - xB) * (yi - xB) / (xi - xB)
    rect_inv = lambda y: (y - xD / (R + 1)) * (R + 1) / R
    strip_inv = lambda y: xB + (y - xB) * (xi - xB) / (yi - xB)
    return rect, strip, rect_inv, strip_inv, xi, yi


def mccabe_thiele(R):
    """Step off ideal stages from the top. Returns (stages incl. reboiler, feed stage, points)."""
    rect, strip, _, _, xi, _ = lines(R)
    pts, x, y, n, feed = [(xD, xD)], xD, xD, 0, None
    while n < 2000:
        xn = x_eq(y)
        n += 1
        pts.append((xn, y))
        if xn <= xB:
            return n - 1 + (x - xB) / (x - xn), feed, pts
        if xn > xi:
            y = rect(xn)
        else:
            feed = feed or n
            y = strip(xn)
        x = xn
        pts.append((x, y))
    raise RuntimeError("stepping did not converge")


R = R_FACTOR * R_min
N_ideal, feed_stage, pts_mt = mccabe_thiele(R)
N_ideal_int = int(np.ceil(N_ideal))
L, V = R * D, (R + 1) * D            # rectifying section, kmol/h
Lb, Vb = L + q * F, V - (1 - q) * F  # stripping section, kmol/h


def gilliland(R):
    """Gilliland correlation (Molokanov equation); returns stages incl. reboiler."""
    X = (R - R_min) / (R + 1)
    Y = 1 - np.exp((1 + 54.4 * X) / (11 + 117.2 * X) * (X - 1) / np.sqrt(X))
    Nm = N_min_fenske + 1
    return (Nm + Y) / (1 - Y)


factors = [1.05, 1.1, 1.2, 1.3, 1.5, 2.0, 3.0]
sens = [(f, f * R_min, mccabe_thiele(f * R_min)[0], gilliland(f * R_min)) for f in factors]

# ----------------------------------------------------------------------------
# 6. Stage efficiency
# ----------------------------------------------------------------------------
T_top, T_bot = T_bubble(xD), T_bubble(xB)
mu_feed = xF * MU_L["EB"] + (1 - xF) * MU_L["SM"]          # cP, at average column temperature
eta_o = 0.492 * (mu_feed * a_mean) ** -0.245               # O'Connell
N_plates_ideal = N_ideal - 1                               # reboiler is one ideal stage
N_actual_oconnell = int(np.ceil(N_plates_ideal / eta_o))


def murphree_trays(R, eta_M):
    """Step off real trays from the reboiler upward on the effective equilibrium curve."""
    rect, strip, rect_inv, strip_inv, xi, yi = lines(R)
    y = y_eq(xB)                                 # reboiler: true equilibrium stage
    pts, n, feed_from_bottom = [(xB, xB), (xB, y)], 0, None
    while n < 2000:
        if y < yi:
            x, y_in = strip_inv(y), y
        else:
            x, y_in = rect_inv(y), y
            feed_from_bottom = feed_from_bottom or n     # trays below the feed tray
        pts.append((x, y_in))
        y_out = y_in + eta_M * (y_eq(x) - y_in)   # y'_e = y + eta_M (y_e - y)
        n += 1
        pts.append((x, y_out))
        if y_out >= xD:
            return n - 1 + (xD - y_in) / (y_out - y_in), feed_from_bottom, pts
        y = y_out
    raise RuntimeError("stepping did not converge")


eta_M = eta_o                                    # design assumption eta_M ~ eta_o
N_trays_murphree, n_below_feed, pts_mur = murphree_trays(R, eta_M)
N_trays_int = int(np.ceil(N_trays_murphree))
feed_tray_actual = N_trays_int - n_below_feed + 1   # numbered from the top


def eta_o_from_murphree(x, L_over_V):
    lam = m_slope(x) / L_over_V                  # m V / L
    return np.log(1 + eta_M * (lam - 1)) / np.log(lam)


eff_check = [
    ("top (x = xD)", xD, L / V),
    ("rectifying, at feed", xF, L / V),
    ("stripping, at feed", xF, Lb / Vb),
    ("bottom (x = xB)", xB, Lb / Vb),
]
eff_rows = [(lab, float(m_slope(x)), float(m_slope(x) / lv), float(eta_o_from_murphree(x, lv))) for lab, x, lv in eff_check]

# ----------------------------------------------------------------------------
# 7. Packed-column alternative: NTU and HTU
# ----------------------------------------------------------------------------
rect, strip, rect_inv, strip_inv, xi, yi = lines(R)
y_reb = y_eq(xB)                                 # vapor leaving the reboiler = vapor entering the packing


def integrand(y, inv):
    return 1.0 / (y_eq(inv(y)) - y)


NOy_rect = quad(integrand, yi, xD, args=(rect_inv,), limit=200)[0]
NOy_strip = quad(integrand, y_reb, yi, args=(strip_inv,), limit=200)[0]
NOy = NOy_rect + NOy_strip

# --- 7.1 flows and physical properties per section ---------------------------
LB_FT2H = 4.8824 / 3600        # 1 lb/ft2.h in kg/m2.s
FT = 0.3048
R_GAS = 8.314


def section_props(x_avg, liq, vap, solute, solvent):
    T = T_bubble(x_avg)
    TK = T + 273.15
    M = x_avg * MW["EB"] + (1 - x_avg) * MW["SM"]
    rho_v = 101325 * M / 1000 / (R_GAS * TK)
    rho_l = x_avg * RHO_L["EB"] + (1 - x_avg) * RHO_L["SM"]
    mu_l = x_avg * MU_L["EB"] + (1 - x_avg) * MU_L["SM"]              # cP
    mu_v = x_avg * MU_V["EB"] + (1 - x_avg) * MU_V["SM"]              # Pa.s
    # Fuller-Schettler-Giddings (vapor), cm2/s, P in atm
    Dv = 1e-3 * TK ** 1.75 * (1 / MW["EB"] + 1 / MW["SM"]) ** 0.5 / (
        1.0 * (FSG_V["EB"] ** (1 / 3) + FSG_V["SM"] ** (1 / 3)) ** 2)
    # Wilke-Chang (liquid), cm2/s: dilute solute in solvent, psi = 1
    Dl = 7.4e-8 * (1.0 * MW[solvent]) ** 0.5 * TK / (MU_L[solvent] * VB[solute] ** 0.6)
    Sc_v = mu_v / (rho_v * Dv * 1e-4)
    Sc_l = mu_l * 1e-3 / (rho_l * Dl * 1e-4)
    return dict(T=T, M=M, rho_v=rho_v, rho_l=rho_l, mu_l=mu_l, mu_v=mu_v, Dv=Dv, Dl=Dl,
                Sc_v=Sc_v, Sc_l=Sc_l, L=liq, V=vap, mL=liq * M / 3600, mV=vap * M / 3600)


sec = {
    "rect": section_props((xF + xD) / 2, L, V, "SM", "EB"),
    "strip": section_props((xF + xB) / 2, Lb, Vb, "EB", "SM"),
}

# --- 7.2 column diameter from the generalized pressure-drop correlation ------
FP = 95.0                   # packing factor, 1.5-in. ceramic Raschig rings (Table 18.1)
Y_FLOOD = 0.10              # Fig. 18.6 ordinate of the top (1.5 in. H2O/ft) curve at X ~ 0.06-0.07
FLOOD_FRACTION = 0.5        # design vapor velocity = 50 % of flooding
G_C = 32.174
for s in sec.values():
    s["X"] = s["mL"] / s["mV"] * np.sqrt(s["rho_v"] / (s["rho_l"] - s["rho_v"]))
    rl, rv = s["rho_l"] / 16.018, s["rho_v"] / 16.018          # lb/ft3
    Gy_flood = np.sqrt(Y_FLOOD * G_C * (rl - rv) * rv / (FP * s["mu_l"] ** 0.1))  # lb/ft2.s
    s["Gy_flood"] = Gy_flood * 3600                             # lb/ft2.h
    s["S_req"] = s["mV"] / (FLOOD_FRACTION * s["Gy_flood"] * LB_FT2H)   # m2
S_col = max(s["S_req"] for s in sec.values())
D_col = np.sqrt(4 * S_col / np.pi)
for s in sec.values():
    s["Gy"] = s["mV"] / S_col / LB_FT2H         # lb/ft2.h
    s["Gx"] = s["mL"] / S_col / LB_FT2H
    s["pct_flood"] = 100 * s["Gy"] / s["Gy_flood"]
    s["Y_design"] = Y_FLOOD * (s["Gy"] / s["Gy_flood"]) ** 2

# --- 7.3 film transfer-unit heights from Figs. 18.21 and 18.22 ----------------
# Fig. 18.22 readings: H_y (ft) for NH3-air-water, 1.5-in. ceramic Raschig rings (below loading)
HY_DATA = np.array([   # Gx, Gy (lb/ft2.h), Hy (ft)
    [500, 200, 1.78], [500, 350, 2.34], [500, 600, 3.02],
    [1500, 200, 0.90], [1500, 350, 1.13], [1500, 600, 1.40],
    [4500, 200, 0.60], [4500, 350, 0.71], [4500, 600, 0.84],
])
A_fit = np.column_stack([np.ones(len(HY_DATA)), np.log(HY_DATA[:, 1]), np.log(HY_DATA[:, 0])])
coef, *_ = np.linalg.lstsq(A_fit, np.log(HY_DATA[:, 2]), rcond=None)
HY_A, HY_B, HY_C = np.exp(coef[0]), coef[1], coef[2]
hy_fit_err = np.max(np.abs(np.exp(A_fit @ coef) / HY_DATA[:, 2] - 1))

SC_NH3 = 0.66               # NH3-air at 25 degC
SC_O2, MU_W = 381.0, 0.894  # O2-water at 25 degC; water viscosity, cP
HX_REF, GX_REF, HX_N = 0.79, 1000.0, 0.2   # Fig. 18.21: H_x = 0.79 ft at Gx = 1000, slope 0.2 (1.5-in.)

for key, inv, (ya, yb) in (("rect", rect_inv, (yi, xD)), ("strip", strip_inv, (y_reb, yi))):
    s = sec[key]
    s["Hy_NH3"] = HY_A * s["Gy"] ** HY_B * s["Gx"] ** HY_C                       # ft
    s["Hy"] = s["Hy_NH3"] * (s["Sc_v"] / SC_NH3) ** 0.5                          # ft
    s["Hx_O2"] = HX_REF * (s["Gx"] / GX_REF) ** HX_N                             # ft
    s["Hx"] = s["Hx_O2"] * (MU_W / s["mu_l"]) ** HX_N * (s["Sc_l"] / SC_O2) ** 0.5   # ft
    s["NOy"] = quad(integrand, ya, yb, args=(inv,), limit=200)[0]
    VL = s["V"] / s["L"]
    # Z = int H_Oy dN_Oy with H_Oy = H_y + m (G_M / L_M) H_x and local slope m
    f = lambda y: (s["Hy"] + m_slope(inv(y)) * VL * s["Hx"]) * integrand(y, inv)
    s["Z"] = quad(f, ya, yb, limit=200)[0] * FT                                   # m
    s["HOy"] = s["Z"] / s["NOy"]                                                  # m (N_Oy-weighted mean)
    s["m_avg"] = (s["HOy"] / FT - s["Hy"]) / (VL * s["Hx"])
Z_total = sec["rect"]["Z"] + sec["strip"]["Z"]
HETP = Z_total / N_plates_ideal

# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
T_feed = T_bubble(xF)

print(f"T_top = {T_top:.2f} C, T_feed = {T_feed:.2f} C, T_bot = {T_bot:.2f} C")
print(f"alpha: D {aD:.4f}, F {aF:.4f}, B {aB:.4f}, mean {a_mean:.4f}")
print(f"D = {D:.2f}, B = {B:.2f} kmol/h")
print(f"y' = {y_p:.4f}, R_min = {R_min:.3f}")
print(f"N_min: Fenske {N_min_fenske:.2f} plates + reboiler; graphical {N_min_graph:.2f} stages")
print(f"R = {R:.3f}: L {L:.1f}, V {V:.1f}, Lbar {Lb:.1f}, Vbar {Vb:.1f} kmol/h; intersection y = {yi:.4f}")
print(f"N_ideal = {N_ideal:.2f} stages incl. reboiler (-> {N_ideal_int}), feed stage {feed_stage}")
for f_, R_, N_, G_ in sens:
    print(f"  R/Rmin {f_:4.2f}  R {R_:6.2f}  N(M-T) {N_:6.1f}  N(Gilliland) {G_:6.1f}")
print(f"mu_feed = {mu_feed:.4f} cP, eta_o (O'Connell) = {eta_o:.4f}, actual trays = {N_actual_oconnell}")
print(f"Murphree stepping: {N_trays_murphree:.2f} trays (-> {N_trays_int}), feed tray {feed_tray_actual} from top")
for row in eff_rows:
    print("  %-22s m = %.4f  mV/L = %.4f  eta_o = %.4f" % row)
print(f"y_reb = {y_reb:.5f}; N_Oy rect {NOy_rect:.2f}, strip {NOy_strip:.2f}, total {NOy:.2f}")
print(f"H_y,NH3 fit: {HY_A:.3f} Gy^{HY_B:.3f} Gx^{HY_C:.3f}, max dev {100*hy_fit_err:.1f} %")
print(f"S = {S_col:.2f} m2, D_col = {D_col:.2f} m")
for k, s in sec.items():
    print(f"[{k}] T {s['T']:.1f} C, rho_v {s['rho_v']:.3f}, rho_l {s['rho_l']:.1f} kg/m3, Sc_v {s['Sc_v']:.3f}, Sc_l {s['Sc_l']:.1f}, "
          f"X {s['X']:.4f}, Gy_flood {s['Gy_flood']:.0f}, Gy {s['Gy']:.0f}, Gx {s['Gx']:.0f} lb/ft2.h")
    print(f"      Hy_NH3 {s['Hy_NH3']:.3f}, Hy {s['Hy']:.3f}, Hx_O2 {s['Hx_O2']:.3f}, Hx {s['Hx']:.3f} ft, m_avg {s['m_avg']:.3f}, "
          f"H_Oy {s['HOy']:.3f} m, N_Oy {s['NOy']:.2f}, Z {s['Z']:.1f} m")
print(f"Z_total = {Z_total:.1f} m, HETP = {HETP:.3f} m")

# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, MUTED, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": "#b5b4ae", "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": "#e6e5e0",
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold", "legend.frameon": False,
    "lines.linewidth": 1.8,
})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=200)
    plt.close(fig)


def stairs(ax, pts, **kw):
    px, py = zip(*pts)
    ax.plot(px, py, color=INK, lw=0.9, label=kw.get("label"), zorder=3)


def base_xy(ax, R=None, legend=True):
    ax.plot(xg, yg, color=BLUE, label="Equilibrium curve")
    ax.plot([0, 1], [0, 1], color=MUTED, lw=1.0, label="y = x")
    if R is not None:
        rc, st, *_ , xi_, yi_ = lines(R)
        ax.plot([xi_, xD], [yi_, xD], color=ORANGE, ls="--", label="Rectifying operating line")
        ax.plot([xB, xi_], [xB, yi_], color=AQUA, ls="-.", label="Stripping operating line")
        ax.plot([xF, xF], [xF, y_eq(xF)], color=YELLOW, ls=":", lw=2.2, label="Feed line (q = 1)")
    ax.set_xlabel("x, mole fraction EB in liquid")
    ax.set_ylabel("y, mole fraction EB in vapor")


# Fig. 1: VLE
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.4))
ax = axs[0]
ax.plot(xg, Tg, color=BLUE, label="Bubble-point curve (T-x)")
ax.plot(yg, Tg, color=ORANGE, ls="--", label="Dew-point curve (T-y)")
ax.set_xlabel("x, y  (mole fraction EB)"); ax.set_ylabel("Temperature, °C")
ax.set_title("(a) T-x-y diagram, EB/SM at 760 mmHg"); ax.legend(loc="upper right")
ax = axs[1]
base_xy(ax)
ax.set_title("(b) x-y diagram (Raoult's law)")
ax.annotate(f"α = {aB:.3f} (x→0) … {aD:.3f} (x→1)", xy=(0.5, y_eq(0.5)), xytext=(0.08, 0.86),
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8), color=INK)
ax.legend(loc="lower right"); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
save(fig, "fig1_vle.png")

# Fig. 2: minimum reflux and total reflux
fig, axs = plt.subplots(1, 2, figsize=(10.5, 5.0))
ax = axs[0]
base_xy(ax)
ax.plot([0, xD], [xD / (R_min + 1), xD], color=ORANGE, ls="--", label=f"Operating line at $R_{{Dm}}$ = {R_min:.2f}")
ax.plot([xF, xF], [xF, y_p], color=YELLOW, ls=":", lw=2.2, label="Feed line (q = 1)")
ax.plot([x_p], [y_p], "o", color=INK, ms=5, zorder=4)
ax.annotate(f"pinch (x', y') = ({x_p:.2f}, {y_p:.4f})", xy=(x_p, y_p), xytext=(0.50, 0.33),
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
ax.annotate(f"intercept $x_D/(R_{{Dm}}+1)$ = {xD/(R_min+1):.4f}", xy=(0, xD / (R_min + 1)), xytext=(0.16, 0.04),
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
ax.set_title("(a) Minimum reflux ratio"); ax.legend(loc="upper left")
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
ax = axs[1]
base_xy(ax)
stairs(ax, pts_total, label=f"Stages at total reflux ({N_min_graph:.1f})")
ax.set_title("(b) Total reflux: minimum number of stages"); ax.legend(loc="lower right")
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
save(fig, "fig2_rmin_nmin.png")

# Fig. 3: McCabe-Thiele, full diagram
fig, ax = plt.subplots(figsize=(7.2, 7.2))
base_xy(ax, R)
stairs(ax, pts_mt, label=f"Ideal stages ({N_ideal_int} incl. reboiler)")
for xv, lab in ((xB, "$x_B$"), (xF, "$x_F$"), (xD, "$x_D$")):
    ax.plot([xv, xv], [0, xv], color=MUTED, lw=0.7, ls=":")
    dx = -0.012 if xv > 0.9 else (0.07 if xv < 0.1 else 0.012)
    ax.text(xv + dx, 0.015, lab, ha="right" if xv > 0.9 else "left", va="bottom")
ax.set_title(f"McCabe-Thiele diagram, $R_D$ = {R_FACTOR}·$R_{{Dm}}$ = {R:.2f}")
ax.legend(loc="upper left"); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
save(fig, "fig3_mccabe_thiele.png")


def number_stages(ax, pts, every, xlim, top_down=True):
    """Label every n-th stage at its corner on the equilibrium curve."""
    corners = pts[1::2]
    for i, (px, py) in enumerate(corners, start=1):
        if i % every == 0 and xlim[0] < px < xlim[1]:
            ax.annotate(str(i), xy=(px, py), xytext=(-3, 3), textcoords="offset points",
                        ha="right", va="bottom", fontsize=8, color=MUTED)


# Fig. 4: McCabe-Thiele, zoomed
fig, axs = plt.subplots(1, 3, figsize=(15, 5.2))
zooms = [("(a) Top of column", (0.88, 1.0)), ("(b) Feed region", (0.33, 0.47)), ("(c) Bottom of column", (0.0, 0.12))]
for ax, (title, (lo, hi)) in zip(axs, zooms):
    base_xy(ax, R)
    stairs(ax, pts_mt, label="Ideal stages")
    number_stages(ax, pts_mt, 5, (lo, hi))
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi + 0.02 if hi < 1 else hi); ax.set_title(title)
axs[1].annotate(f"feed stage {feed_stage}", xy=(xF, yi), xytext=(0.40, 0.385),
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
axs[0].legend(loc="upper left", fontsize=8.5)
save(fig, "fig4_mccabe_thiele_zoom.png")

# Fig. 5: reflux sensitivity
ff = np.linspace(1.04, 3.0, 60)
fig, ax = plt.subplots(figsize=(7.2, 4.4))
ax.plot(ff, [mccabe_thiele(f * R_min)[0] for f in ff], color=BLUE, label="McCabe-Thiele stepping")
ax.plot(ff, [gilliland(f * R_min) for f in ff], color=ORANGE, ls="--", label="Gilliland correlation (Molokanov eq.)")
ax.axhline(N_min_graph, color=MUTED, lw=0.9, ls=":")
ax.text(2.98, N_min_graph + 2, f"$N_{{min}}$ = {N_min_graph:.1f} stages", ha="right", color=MUTED)
ax.plot([R_FACTOR], [N_ideal], "o", color=INK, ms=6, zorder=4)
ax.annotate(f"design: {R_FACTOR}·$R_{{Dm}}$, {N_ideal:.1f} stages", xy=(R_FACTOR, N_ideal), xytext=(1.45, N_ideal + 18),
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
ax.set_xlabel("$R_D / R_{Dm}$"); ax.set_ylabel("Ideal stages (incl. reboiler)")
ax.set_title("Number of ideal stages vs. reflux ratio"); ax.legend(loc="upper right")
ax.set_xlim(1.0, 3.0); ax.set_ylim(0, None)
save(fig, "fig5_reflux_sensitivity.png")

# Fig. 6: Murphree efficiency, effective equilibrium curve
def y_eff(x):
    yop = np.where(x >= xF, rect(x), strip(x))
    return yop + eta_M * (ye_spline(x) - yop)


xs = np.linspace(xB, xD, 600)
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.4))
for ax, (lo, hi), title in zip(axs, ((0.0, 1.0), (0.33, 0.47)), ("(a) Full diagram", "(b) Feed region")):
    base_xy(ax, R)
    ax.plot(xs, y_eff(xs), color=INK, ls=(0, (4, 2)), lw=1.3, label=f"Effective equilibrium curve ($\\eta_M$ = {eta_M:.2f})")
    stairs(ax, pts_mur, label=f"Actual trays ({N_trays_int}) + reboiler")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi + (0.02 if hi < 1 else 0)); ax.set_title(title)
axs[0].set_aspect("equal"); axs[0].legend(loc="upper left", fontsize=8.5)
save(fig, "fig6_murphree.png")

# Fig. 7: NTU integrand
fig, ax = plt.subplots(figsize=(7.6, 4.4))
ys = np.linspace(y_reb, yi, 300); yr = np.linspace(yi, xD, 300)
fs = [integrand(y, strip_inv) for y in ys]; fr = [integrand(y, rect_inv) for y in yr]
ax.plot(ys, fs, color=AQUA); ax.fill_between(ys, fs, color=AQUA, alpha=0.25, label=f"Stripping section: $N_{{Oy}}$ = {NOy_strip:.1f}")
ax.plot(yr, fr, color=ORANGE); ax.fill_between(yr, fr, color=ORANGE, alpha=0.25, label=f"Rectifying section: $N_{{Oy}}$ = {NOy_rect:.1f}")
ax.axvline(yi, color=MUTED, lw=0.8, ls=":"); ax.text(yi, ax.get_ylim()[1] * 0.97, " feed", color=MUTED, va="top")
ax.set_xlabel("y, mole fraction EB in vapor"); ax.set_ylabel("1 / (y* − y)")
ax.set_title("Number of transfer units: $N_{Oy}$ = ∫ dy / (y* − y)"); ax.legend(loc=(0.48, 0.45))
ax.set_ylim(0, None)
save(fig, "fig7_ntu.png")
print("figures written to", OUT)
