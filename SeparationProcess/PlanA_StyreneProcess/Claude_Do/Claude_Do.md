# Distillation Column Design: Ethylbenzene / Styrene at 1 atm
 
Separation Process term project. A continuous distillation column separating ethylbenzene (EB, light component) from styrene monomer (SM) is designed by hand-calculation methods from the lecture notes, implemented in Python ([`ebsm_design.py`](ebsm_design.py)).
 
Lecture references are written as **L6 p.26** = `SeparationProcess6.pdf`, slide 26.
 
The calculation runs in dependency order: VLE, material balance, minimum reflux, minimum stages, McCabe–Thiele, stage efficiency, and finally NTU/HTU. The transfer-unit calculation comes last because the NTU integral needs the operating lines, which exist only after the reflux ratio is fixed.
 
## Summary of results
 
| Quantity | Result | Method |
|---|---|---|
| Relative volatility (top / feed / bottom) | 1.278 / 1.271 / 1.266, mean 1.272 | Antoine + Raoult |
| Distillate / bottoms | 39.8 / 60.2 kmol/h per 100 kmol/h feed | Material balance |
| Minimum reflux ratio $R_{Dm}$ | 9.05 | Pinch at the feed line |
| Minimum number of stages $N_{min}$ | 37.2 plates + reboiler | Fenske |
| Operating reflux ratio $R_D$ | 10.86 ($1.2\,R_{Dm}$) | Design choice |
| Ideal stages | 78 (77 plates + reboiler), feed on stage 38 | McCabe–Thiele |
| Overall efficiency $\eta_o$ | 0.66 | O'Connell |
| Actual trays | 115 + reboiler, feed on tray 57 | Murphree effective equilibrium curve |
| Transfer units $N_{Oy}$ | 76.3 (38.1 rectifying + 38.2 stripping) | $\int dy/(y^{\ast}-y)$ |
| Transfer-unit height $H_{Oy}$ | 0.78 m (rectifying), 0.72 m (stripping) | Figs. 18.21, 18.22 |
| Packed height | 57 m of 1½-in. ceramic Raschig rings (HETP ≈ 0.75 m) | $Z = H_{Oy} N_{Oy}$ |
 
## 1. Design basis
 
| Item | Value | Status |
|---|---|---|
| System | Ethylbenzene (A, light) / styrene (B, heavy) | Given |
| Pressure | 1 atm (760 mmHg), uniform in the column | Given |
| Polymerization of styrene | Ignored | Given |
| Feed condition | Saturated liquid, $q = 1$ | Given |
| Distillate composition $x_D$ | 0.99 mole fraction EB | Given |
| Feed composition $x_F$ | 0.40 mole fraction EB | Assumed |
| Bottoms composition $x_B$ | 0.01 mole fraction EB (99 % styrene) | Assumed |
| Feed rate $F$ | 100 kmol/h (calculation basis) | Assumed |
| Reflux ratio | $R_D = 1.2\,R_{Dm}$ (as in Example 21.5, L7 p.19) | Assumed |
| Condenser / reboiler | Total condenser; partial reboiler = one ideal stage | Assumed |
 
All assumed values are constants at the top of the script.
 
Assumptions of the McCabe–Thiele method (L6 p.12): equal molar latent heats (constant molal overflow), negligible sensible heat and heat of mixing, adiabatic column, uniform pressure. The molar latent heats at the normal boiling points are 35.6 kJ/mol (EB) and 37.1 kJ/mol (SM), within about 4 % of each other, so constant molal overflow is a good approximation for this system.
 
## 2. Vapor–liquid equilibrium
 
Vapor pressures come from the Antoine equation, $\log_{10} P^{sat}\,[\mathrm{mmHg}] = A - B/(T\,[^{\circ}\mathrm{C}] + C)$:
 
| Component | A | B | C | Normal boiling point |
|---|---|---|---|---|
| Ethylbenzene | 6.95719 | 1424.255 | 213.206 | 136.2 °C |
| Styrene | 6.92409 | 1420.0 | 206.0 | 145.2 °C |
 
Both components are C8 aromatics, so the ideal-gas, ideal-solution model is used (L1 p.33, L8 p.2):
 
```math
K_i = \frac{y_i}{x_i} = \frac{P_i^{sat}}{P}, \qquad
\alpha_{AB} = \frac{K_A}{K_B} = \frac{P_A^{sat}}{P_B^{sat}}
```
 
For each liquid composition $x$, the bubble-point temperature is found from $\sum K_i x_i = 1$ (L8 p.3–4), which gives $y = K_A x$. The dew-point calculation $\sum y_i/K_i = 1$ gives $x$ from $y$. Repeating over $0 \le x \le 1$ gives the T-x-y and x-y diagrams.
 
![T-x-y and x-y diagrams](images/fig1_vle.png)
 
| Location | $x$ | Bubble point, °C | $\alpha$ |
|---|---|---|---|
| Distillate | 0.99 | 136.27 | 1.2778 |
| Feed | 0.40 | 141.36 | 1.2711 |
| Bottoms | 0.01 | 145.10 | 1.2664 |
 
The mean relative volatility is the geometric mean of the three values (L8 p.12):
 
```math
\bar{\alpha} = \sqrt[3]{\alpha_D\,\alpha_F\,\alpha_B} = 1.2718
```
 
The boiling points differ by only 9 °C and the equilibrium curve lies close to the diagonal, so this is a difficult separation. All graphical constructions below use the rigorous Raoult's-law curve; $\bar{\alpha}$ is used only in the Fenske and O'Connell equations.
 
## 3. Material balance
 
Overall and component balances (L6 p.9):
 
```math
F = D + B, \qquad F x_F = D x_D + B x_B
```
 
```math
\frac{D}{F} = \frac{x_F - x_B}{x_D - x_B} = \frac{0.40 - 0.01}{0.99 - 0.01} = 0.398
```
 
$D = 39.80$ kmol/h and $B = 60.20$ kmol/h.
 
## 4. Minimum reflux ratio
 
At minimum reflux the rectifying operating line passes through $(x_D, x_D)$ and the point $(x', y')$ where the feed line meets the equilibrium curve (L6 p.27). For $q = 1$ the feed line is vertical, so $x' = x_F = 0.40$ and $y' = y_e(0.40) = 0.4587$.
 
```math
R_{Dm} = \frac{x_D - y'}{y' - x'} = \frac{0.99 - 0.4587}{0.4587 - 0.40} = 9.05
```
 
The equilibrium curve has no inflection, so the pinch is at the feed and this equation applies (L6 p.28, case a).
 
## 5. Minimum number of stages
 
At total reflux both operating lines coincide with the diagonal. The Fenske equation (L6 p.26) gives
 
```math
N_{min} = \frac{\ln\left[\dfrac{x_D (1 - x_B)}{x_B (1 - x_D)}\right]}{\ln \bar{\alpha}} - 1
= \frac{\ln\left[\dfrac{(0.99)(0.99)}{(0.01)(0.01)}\right]}{\ln 1.2718} - 1
= \frac{9.190}{0.2404} - 1 = 37.2
```
 
That is 37.2 ideal plates plus the reboiler (38.2 stages). Stepping between the diagonal and the rigorous equilibrium curve gives the same 38.2 stages.
 
![Minimum reflux ratio and total reflux](images/fig2_rmin_nmin.png)
 
## 6. McCabe–Thiele diagram
 
### 6.1 Operating lines
 
Reflux ratio: $R_D = 1.2 \times 9.05 = 10.86$.
 
Internal flows for constant molal overflow (L6 p.14, p.21):
 
| Stream | Equation | kmol/h |
|---|---|---|
| $L$ (rectifying liquid) | $R_D D$ | 432.2 |
| $V$ (rectifying vapor) | $(R_D + 1) D$ | 472.0 |
| $\bar{L}$ (stripping liquid) | $L + qF$ | 532.2 |
| $\bar{V}$ (stripping vapor) | $V - (1-q)F$ | 472.0 |
 
Rectifying operating line (L6 p.14):
 
```math
y = \frac{R_D}{R_D + 1}\,x + \frac{x_D}{R_D + 1} = 0.9157\,x + 0.0835
```
 
Feed line (L6 p.21): $q = 1$, so it is the vertical line $x = x_F = 0.40$. It meets the rectifying line at $(0.40,\ 0.4497)$.
 
Stripping operating line (L6 p.18), through $(x_B, x_B)$ and the intersection point:
 
```math
y = \frac{\bar{L}}{\bar{V}}\,x - \frac{B x_B}{\bar{V}} = 1.1275\,x - 0.00128
```
 
### 6.2 Number of ideal stages
 
Stages are stepped off from $(x_D, x_D)$ downward. The construction switches to the stripping line at the first stage whose liquid composition falls below the intersection, which is the optimum feed location (L6 p.23).
 
![McCabe-Thiele diagram](images/fig3_mccabe_thiele.png)
 
![McCabe-Thiele diagram, enlarged at the top, feed and bottom](images/fig4_mccabe_thiele_zoom.png)
 
Result: 77.2 stages, rounded up to **78 ideal stages = 77 ideal plates + reboiler**, with the feed on **stage 38** from the top.
 
### 6.3 Effect of reflux ratio
 
| $R_D/R_{Dm}$ | $R_D$ | Ideal stages (McCabe–Thiele) | Ideal stages (Gilliland) |
|---|---|---|---|
| 1.05 | 9.50 | 102.0 | 101.5 |
| 1.10 | 9.96 | 89.5 | 90.6 |
| **1.20** | **10.86** | **77.2** | **77.9** |
| 1.30 | 11.77 | 70.3 | 70.3 |
| 1.50 | 13.58 | 62.3 | 61.6 |
| 2.00 | 18.10 | 53.4 | 52.4 |
| 3.00 | 27.15 | 47.0 | 46.2 |
 
The Gilliland correlation (Molokanov equation, L8 p.24) needs only $N_{min}$ and $R_{Dm}$ and agrees with the graphical result within 2 %, which checks sections 4 to 6 against each other.
 
![Number of ideal stages versus reflux ratio](images/fig5_reflux_sensitivity.png)
 
The number of stages rises steeply below about $1.2\,R_{Dm}$, while the vapor rate, column diameter and heat duties grow in proportion to $R_D + 1$ (L6 p.29). $1.2\,R_{Dm}$ sits at the knee of the curve.
 
## 7. Stage efficiency
 
### 7.1 Overall efficiency (O'Connell)
 
The O'Connell correlation (L7 p.15) uses the liquid viscosity of the feed at the average column temperature, here 141 °C:
 
```math
\mu_L = 0.40\,(0.227) + 0.60\,(0.234) = 0.231\ \mathrm{cP}
```
 
```math
\eta_o = 0.492\,(\mu_L \bar{\alpha})^{-0.245} = 0.492\,(0.231 \times 1.272)^{-0.245} = 0.664
```
 
Number of actual trays (L7 p.9), with the reboiler counted as one ideal stage:
 
```math
N_{actual} = \frac{N_{ideal}}{\eta_o} = \frac{77.2 - 1}{0.664} = 114.7 \rightarrow 115\ \text{trays}
```
 
### 7.2 Murphree efficiency on the McCabe–Thiele diagram
 
For a column with both rectifying and stripping sections, $\eta_M \approx \eta_o$ can be assumed in design (L7 p.14), so $\eta_M = 0.664$. The true equilibrium curve is replaced by the effective equilibrium curve (L7 p.12)
 
```math
y_e' = y + \eta_M\,(y_e - y)
```
 
where $y$ is the operating-line value at the same $x$. The reboiler keeps the true equilibrium curve. Stepping from the reboiler upward gives 114.8, rounded up to **115 actual trays plus the reboiler**, with the feed on **tray 57** from the top.
 
![Effective equilibrium curve and actual trays](images/fig6_murphree.png)
 
### 7.3 Consistency of the two efficiencies
 
The relation between the two efficiencies for straight equilibrium and operating lines (L7 p.14) is
 
```math
\eta_o = \frac{\ln\left[1 + \eta_M\,(mV/L - 1)\right]}{\ln(mV/L)}
```
 
with $m$ the local slope of the equilibrium curve.
 
| Location | $m$ | $mV/L$ | $\eta_o$ from $\eta_M = 0.664$ |
|---|---|---|---|
| Top ($x = x_D$) | 0.786 | 0.858 | 0.647 |
| Rectifying section at the feed | 1.037 | 1.132 | 0.678 |
| Stripping section at the feed | 1.037 | 0.920 | 0.655 |
| Bottom ($x = x_B$) | 1.260 | 1.117 | 0.676 |
 
$mV/L$ stays between 0.86 and 1.13, so $\eta_o$ stays within 3 % of $\eta_M$. This is why both routes give 115 trays.
 
## 8. Packed-column alternative: NTU and HTU
 
The same separation is sized as a packed column with $Z = H_{Oy} N_{Oy}$ (L3 p.19). The packing is 1½-in. ceramic Raschig rings, the packing for which the lecture notes give both film correlations (L3 p.27–28).
 
### 8.1 Number of transfer units
 
```math
N_{Oy} = \int_{y_b}^{y_a} \frac{dy}{y^{\ast} - y}
```
 
Here $y$ follows the operating line and $y^{\ast} = y_e(x)$ is the vapor in equilibrium with the liquid at the same height. With constant molal overflow the vapor rate is constant, so this form holds without a dilute-gas correction. The integral is evaluated numerically for each section:
 
| Section | Limits in $y$ | $N_{Oy}$ |
|---|---|---|
| Stripping | 0.0126 (vapor from the reboiler) to 0.4497 (feed) | 38.2 |
| Rectifying | 0.4497 (feed) to 0.99 (top) | 38.1 |
| Total | | **76.3** |
 
![Integrand of the NTU integral](images/fig7_ntu.png)
 
The driving force $y^{\ast} - y$ becomes small at both ends of the column, where most of the transfer units are spent. $N_{Oy}$ is nearly equal to the 76.2 ideal plates because the operating and equilibrium lines are almost parallel ($mV/L \approx 1$; L3 p.19).
 
### 8.2 Flow rates and physical properties
 
Properties are evaluated at the mean composition of each section. Pure-component values at 140 °C are taken from the `thermo` Python library.
 
| Property | Rectifying | Stripping |
|---|---|---|
| Mean temperature, °C | 138.7 | 143.2 |
| Vapor mass flow, kg/s | 13.84 | 13.71 |
| Liquid mass flow, kg/s | 12.67 | 15.46 |
| Vapor density $\rho_y$, kg/m³ (ideal gas) | 3.12 | 3.06 |
| Liquid density $\rho_x$, kg/m³ | 767 | 783 |
| Liquid viscosity, cP | 0.229 | 0.233 |
| Vapor viscosity, μPa·s | 8.88 | 8.84 |
| Vapor diffusivity $D_v$, cm²/s | 0.0507 | 0.0517 |
| Liquid diffusivity $D_v$, cm²/s | 7.4 × 10⁻⁵ | 6.9 × 10⁻⁵ |
| $Sc$ (vapor) | 0.561 | 0.559 |
| $Sc$ (liquid) | 40.6 | 42.8 |
 
Diffusivities use the correlations from Lecture 2. Vapor: Fuller–Schettler–Giddings (L2, slide 17) with diffusion volumes 131.6 (EB) and 127.6 (SM). Liquid: Wilke–Chang (L2, slide 19) with $\psi = 1$ and molar volumes at the boiling point of 139.6 (EB) and 132.8 (SM) cm³/mol; the minor component of each section is taken as the solute.
 
### 8.3 Column cross-section
 
The mass velocities needed for the HTU correlations follow from the generalized pressure-drop correlation, Fig. 18.6 (L3 p.13).
 
| Step | Rectifying | Stripping |
|---|---|---|
| Flow parameter $\dfrac{G_x}{G_y}\sqrt{\dfrac{\rho_y}{\rho_x - \rho_y}}$ | 0.059 | 0.071 |
| Packing factor $F_p$ (Table 18.1, L3 p.8) | 95 | 95 |
| Ordinate at flooding (chart reading) | 0.10 | 0.10 |
| $G_y$ at flooding, lb/ft²·h | 2174 | 2174 |
| Design $G_y$, lb/ft²·h | 1087 | 1077 |
| Design $G_x$, lb/ft²·h | 995 | 1214 |
 
The flooding pressure drop (about 2 in. H₂O/ft, L3 p.11) lies above the highest curve in Fig. 18.6, so the 1.5 in. H₂O/ft curve is used as a conservative flooding limit. The design vapor velocity is half the flooding velocity (L3 p.12), which gives
 
```math
S = \frac{13.84\ \mathrm{kg/s}}{0.5 \times 2.95\ \mathrm{kg/m^2 s}} = 9.39\ \mathrm{m^2}, \qquad D_c = 3.46\ \mathrm{m}
```
 
At the design point the chart ordinate is 0.025, which corresponds to a pressure drop of about 0.2 in. H₂O/ft.
 
### 8.4 Height of a transfer unit
 
**Gas film.** Fig. 18.22 (L3 p.28) gives $H_y$ for ammonia absorption in water on 1½-in. ceramic Raschig rings. Nine points read from the dashed $H_y$ curves below the loading point:
 
| $G_x$ (rows), $G_y$ (columns), lb/ft²·h | 200 | 350 | 600 |
|---|---|---|---|
| 500 | 1.78 | 2.34 | 3.02 |
| 1500 | 0.90 | 1.13 | 1.40 |
| 4500 | 0.60 | 0.71 | 0.84 |
 
A least-squares fit to these readings (maximum deviation 10 %), followed by the Schmidt-number correction from the same slide:
 
```math
H_{y,\mathrm{NH_3}} = 6.27\,G_y^{0.397}\,G_x^{-0.540}\ \mathrm{ft}, \qquad
H_y = H_{y,\mathrm{NH_3}} \left(\frac{Sc}{0.66}\right)^{1/2}
```
 
The fitted exponents fall in the ranges stated on the slide (0.3 to 0.4 for $G_y$; −0.7 to −0.4 for $G_x$).
 
**Liquid film.** Fig. 18.21 (L3 p.27) gives $H_x$ for oxygen desorption from water at 25 °C. For 1½-in. rings the line reads $H_x = 0.79$ ft at $G_x = 1000$ lb/ft²·h with a slope of 0.2. The slide's equation $H_x = \frac{1}{\alpha}(G_x/\mu)^n Sc^{0.5}$ is applied as a ratio to that reference system (water at 25 °C: $\mu$ = 0.894 cP; $Sc$ = 381 for O₂ in water, ref. 2, Ch. 18):
 
```math
H_x = 0.79 \left(\frac{G_x}{1000}\right)^{0.2} \left(\frac{0.894}{\mu_L}\right)^{0.2} \left(\frac{Sc}{381}\right)^{1/2}\ \mathrm{ft}
```
 
**Overall.** The film heights are combined with the local slope $m$ of the equilibrium curve (L3 p.21):
 
```math
H_{Oy} = H_y + m\,\frac{G_M}{L_M}\,H_x
```
 
| Quantity | Rectifying | Stripping |
|---|---|---|
| $H_{y,\mathrm{NH_3}}$, ft | 2.41 | 2.16 |
| $H_y$, ft | 2.22 | 1.99 |
| $H_x$ (O₂–water), ft | 0.79 | 0.82 |
| $H_x$, ft | 0.34 | 0.36 |
| Mean slope $m$ | 0.88 | 1.16 |
| $G_M/L_M = V/L$ | 1.092 | 0.887 |
| $H_{Oy}$, ft | 2.55 | 2.36 |
| $H_{Oy}$, m | **0.78** | **0.72** |
 
About 85 % of the resistance is in the vapor film.
 
### 8.5 Packed height
 
Because $m$ varies along the column, the height is integrated as $Z = \int H_{Oy}\,dN_{Oy}$; the $H_{Oy}$ values above are the resulting section averages.
 
| Section | $N_{Oy}$ | $H_{Oy}$, m | $Z$, m |
|---|---|---|---|
| Rectifying | 38.1 | 0.78 | 29.6 |
| Stripping | 38.2 | 0.72 | 27.4 |
| Total | 76.3 | | **57.1** |
 
The equivalent HETP is $57.1 / 76.2 = 0.75$ m.
 
Table 18.1 (L3 p.8) lists the relative mass-transfer coefficient $f_p$ of other packings against 1½-in. Raschig rings ($f_p = 1.0$). For 1½-in. metal Pall rings, $f_p = 1.36$ and $F_p = 40$: at the same mass velocities $H_{Oy}$ drops to about 0.55 m and the packed height to about 42 m.
 
## 9. Discussion and limitations
 
- **The column is very tall.** 115 trays or 57 m of packing is the direct result of $\bar{\alpha} = 1.27$ combined with 99 % purity at both ends. With the same VLE model, $\alpha$ rises to about 1.44 at 50 mmHg, so vacuum operation would shorten the column considerably, in addition to limiting polymerization.
- **Uniform pressure is not realistic for the tray column.** At the typical 50–100 mm H₂O per tray (L7 p.6), 115 trays give a pressure drop of 0.6–1.1 bar, which would raise the bottom temperature and lower $\alpha$ there. The packed column loses about 9 kPa (0.2 in. H₂O/ft over 57 m). The stage counts above neglect both effects.
- **The HTU values are order-of-magnitude estimates.** The correlations are for aqueous systems on ceramic Raschig rings, the chart readings carry roughly ±10 % uncertainty, and the design $G_y$ (about 1090 lb/ft²·h) is slightly beyond the plotted range of Fig. 18.22 (1000 lb/ft²·h). The fit uses only points below loading, which overestimates $H_y$ near loading and is therefore conservative.
- **Antoine constants** are common handbook values. They agree with the vapor-pressure correlations in `thermo` within 0.6 % between 100 and 145 °C, but the source should be replaced by the one cited in the final report (for example Perry or the NIST WebBook).
- **Ideal solution.** Raoult's law is used without activity coefficients; a check against a simulator (COCO or DWSIM) is the natural next step.
- **Assumed specifications.** $x_F$, $x_B$, $F$ and the factor $R_D/R_{Dm}$ were not specified and should be confirmed.
## How to reproduce
 
```bash
pip install numpy scipy matplotlib
python ebsm_design.py
```
 
The script prints every number quoted above and writes the seven figures to `images/`.
 
```
.
├── README.md
├── ebsm_design.py
└── images/
    ├── fig1_vle.png
    ├── fig2_rmin_nmin.png
    ├── fig3_mccabe_thiele.png
    ├── fig4_mccabe_thiele_zoom.png
    ├── fig5_reflux_sensitivity.png
    ├── fig6_murphree.png
    └── fig7_ntu.png
```
 
## References
 
1. Lecture notes, Separation Process (SeparationProcess1–9.pdf, SeparationProcess_ColumnDesign.pdf).
2. W. L. McCabe, J. C. Smith, P. Harriott, *Unit Operations of Chemical Engineering*, 7th ed., McGraw-Hill, 2005. Chapters 17, 18, 21, 22.
3. J. D. Seader, E. J. Henley, D. K. Roper, *Separation Process Principles*, 3rd ed., Wiley.
4. C. Bell and contributors, `thermo`: Chemical properties component of Chemical Engineering Design Library (ChEDL), used for density, viscosity and latent-heat values.
 
