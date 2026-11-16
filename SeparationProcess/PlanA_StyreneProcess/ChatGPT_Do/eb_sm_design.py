#!/usr/bin/env python3
"""Atmospheric EB/SM preliminary distillation design.

Run: python eb_sm_design.py
Dependencies: numpy, scipy, matplotlib. No network or external property package.
All compositions are ethylbenzene mole fractions. F, L, V, D, B: kmol/h.
Raoult VLE, ideal vapor, constant molar overflow, q=1, total condenser,
equilibrium partial reboiler, no polymerization, uniform 101325 Pa.

Property coefficients were extracted from chemicals 1.5.2:
  vapor_pressure.Psat_data_Perrys2_8 (Perry 8e, Table 2-8; DIPPR 101);
  viscosity.mu_data_VDI_PPDS_7 (VDI Heat Atlas 2010; PPDS9).
The equations are independently implemented below; chemicals is not required.
HOG is an explicit scenario input, NOT an experimentally predicted packing HTU.
O'Connell overall efficiency and Murphree efficiency are distinct quantities;
the Murphree calculation uses their approximate equality as a separate assumption.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from functools import lru_cache
from math import ceil, exp, log, sqrt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# ln(Psat / Pa) = A + B/T + C*ln(T) + D*T**E; T in K.
PSAT = {
    "EB": (89.063, -7733.7, -9.917, 5.986e-6, 2.0),
    "SM": (105.93, -8685.9, -12.42, 7.5583e-6, 2.0),
}
PSAT_RANGE = {"EB": (178.2, 617.15), "SM": (242.54, 636.0)}
# PPDS9: mu/Pa.s = E*exp(A*X**(1/3)+B*X**(4/3)); X=(C-T)/(T-D).
MU = {
    "EB": (1.62469, 0.79236, 745.185, -16.276, 2.899e-5),
    "SM": (2.52817, 0.19295, 1189.84, 44.635, 5.38e-6),
}


@dataclass(frozen=True)
class Basis:
    F: float = 100.0
    zF: float = 0.50
    xD: float = 0.99
    xB: float = 0.01  # provisional SM purity specification, 99 mol% SM
    P: float = 101325.0
    reflux_factor: float = 1.50
    HOG: float = 0.50  # assumed packing HTU, metres


def psat(T, component):
    A, B, C, D, E = PSAT[component]
    return exp(A + B / T + C * log(T) + D * T**E)


def viscosity(T, component):
    A, B, C, D, E = MU[component]
    X = (C - T) / (T - D)
    root = float(np.cbrt(X))
    return E * exp(A * root + B * X * root)  # Pa.s


@lru_cache(maxsize=60000)
def bubble_T(x, P):
    if not 0.0 <= x <= 1.0:
        raise ValueError("Mole fraction outside [0,1]")
    T = brentq(lambda t: x * psat(t, "EB") + (1 - x) * psat(t, "SM") - P,
               300.0, 500.0, xtol=2e-11)
    for lo, hi in PSAT_RANGE.values():
        if not lo <= T <= hi:
            raise ValueError("Vapor pressure correlation extrapolated")
    return T


def equilibrium_y(x, P):
    return float(x * psat(bubble_T(float(x), P), "EB") / P)


def equilibrium_x(y, P):
    return brentq(lambda x: equilibrium_y(x, P) - y, 0.0, 1.0, xtol=5e-13)


def relative_volatility(x, P):
    T = bubble_T(float(x), P)
    return psat(T, "EB") / psat(T, "SM")


def operating_data(basis, reflux_factor=None):
    b = basis
    if not (0 < b.xB < b.zF < b.xD < 1 and b.F > 0 and b.HOG > 0):
        raise ValueError("Require 0<xB<zF<xD<1, F>0, HOG>0")
    D = b.F * (b.zF - b.xB) / (b.xD - b.xB)
    B = b.F - D
    yp = equilibrium_y(b.zF, b.P)
    Rmin = (b.xD - yp) / (yp - b.zF)
    factor = b.reflux_factor if reflux_factor is None else reflux_factor
    if factor <= 1:
        raise ValueError("Finite-stage design requires reflux_factor > 1")
    R = factor * Rmin
    L, V = R * D, (R + 1) * D
    Ls, Vs = L + b.F, V  # q=1
    ar, br, ass, bs = L / V, D * b.xD / V, Ls / Vs, -B * b.xB / Vs
    return dict(D=D, B=B, Rmin=Rmin, R=R, L=L, V=V, Ls=Ls, Vs=Vs,
                ar=ar, br=br, ass=ass, bs=bs, yp=yp, yI=ar*b.zF+br)


def op_y(x, basis, d):
    return d["ar"]*x+d["br"] if x >= basis.zF else d["ass"]*x+d["bs"]


def equilibrium_steps(basis, d=None, total_reflux=False, eq_override=None):
    """Top-down McCabe steps; reboiler is INCLUDED in total stage count.

    The last fractional stage uses linear interpolation in liquid composition.
    Integer rounding gives a preliminary capacity count, not an exact solution
    of all specifications for a fixed integer-stage MESH model.
    """
    b = basis
    eq = eq_override or (lambda x: equilibrium_y(x, b.P))
    y, old_x = b.xD, b.xD
    points, rows, feed = [(b.xD, b.xD)], [], None
    for n in range(1, 20001):
        x = brentq(lambda z: eq(z) - y, 0., 1., xtol=5e-13)
        if feed is None and x <= b.zF:
            feed = n
        terminal = x <= b.xB
        next_y = x if total_reflux else op_y(x, b, d)
        points.append((x, y))
        rows.append(dict(stage=n, x=x, y=y, y_next=next_y,
                         section="rectifying" if x > b.zF else "stripping",
                         terminal_step=terminal))
        if terminal:
            fraction = (old_x-b.xB)/(old_x-x)
            return dict(N_fractional=n-1+fraction, N_integer=n, feed_stage=feed,
                        rows=rows, points=points)
        if x >= old_x - 1e-13:
            raise RuntimeError("Pinch or infeasible operating line")
        points.append((x, next_y))
        old_x, y = x, next_y
    raise RuntimeError("Stage iteration did not converge")


def packed_ntu(basis, d):
    """Overall gas-phase NTU for packing only, excluding partial reboiler.

    Bottom vapor is yR*=yeq(xB), not xB. Liquid leaving the packing is
    x_pack,out=(V*yR*+B*xB)/L_s. The feed plane is the operating-line
    intersection. The two integrals use y*(x)-y > 0 for distillation.
    """
    b = basis
    y_reboiler = equilibrium_y(b.xB, b.P)
    x_pack_bottom = (y_reboiler - d["bs"]) / d["ass"]
    r = lambda x: d["ar"] * x + d["br"]
    s = lambda x: d["ass"] * x + d["bs"]
    Nr, er = quad(lambda x: d["ar"]/(equilibrium_y(x,b.P)-r(x)),
                  b.zF, b.xD, epsabs=1e-9, epsrel=1e-10)
    Ns, es = quad(lambda x: d["ass"]/(equilibrium_y(x,b.P)-s(x)),
                  x_pack_bottom, b.zF, epsabs=1e-9, epsrel=1e-10)
    # Independent change-of-variable check, directly integrating in y.
    Nr_y = quad(lambda y: 1/(equilibrium_y((y-d["br"])/d["ar"], b.P)-y),
                d["yI"], b.xD, epsabs=1e-9)[0]
    Ns_y = quad(lambda y: 1/(equilibrium_y((y-d["bs"])/d["ass"], b.P)-y),
                y_reboiler, d["yI"], epsabs=1e-9)[0]
    return dict(NOG_R=Nr, NOG_S=Ns, NOG_total=Nr+Ns,
                y_reboiler=y_reboiler, x_pack_bottom=x_pack_bottom,
                Z_R=b.HOG*Nr, Z_S=b.HOG*Ns, Z_total=b.HOG*(Nr+Ns),
                integration_error=max(er,es),
                integration_variable_check=max(abs(Nr-Nr_y),abs(Ns-Ns_y)))


def murphree_steps(basis, d, eta_M):
    """A separate approximate case with eta_M set equal to estimated eta_o.

    Uses an effective equilibrium curve only for trays, then a TRUE equilibrium
    reboiler. Thus no tray efficiency is applied to the reboiler.
    """
    b = basis
    effective = lambda x: eta_M*equilibrium_y(x,b.P)+(1-eta_M)*op_y(x,b,d)
    y_reboiler = equilibrium_y(b.xB,b.P)
    y, previous_inlet, feed = b.xD, b.xD, None
    points, rows = [(b.xD,b.xD)], []
    for n in range(1,30001):
        x = brentq(lambda xx: effective(xx)-y, 0., 1., xtol=5e-13)
        yi = op_y(x,b,d)
        if feed is None and x <= b.zF:
            feed = n
        points.extend([(x,y),(x,yi)])
        rows.append(dict(tray=n,x=x,y_out=y,y_in=yi,
                         eta_M=(y-yi)/(equilibrium_y(x,b.P)-yi)))
        if yi <= y_reboiler:
            fraction=(previous_inlet-y_reboiler)/(previous_inlet-yi)
            # Endpoint overshoot is only for graphical rounding, not a solved
            # product composition at fixed D, R, and an integer tray count.
            points.append((equilibrium_x(yi,b.P),yi))
            return dict(N_trays=n, N_trays_fractional=n-1+fraction,
                        feed_tray=feed, eta_M=eta_M, rows=rows, points=points)
        previous_inlet, y = yi, yi
    raise RuntimeError("Murphree stepping did not converge")


def design(basis):
    b=basis
    d=operating_data(b)
    stages=equilibrium_steps(b,d)
    minimum=equilibrium_steps(b,total_reflux=True)
    ntu=packed_ntu(b,d)
    TD, TF, TB = (bubble_T(x,b.P) for x in [b.xD,b.zF,b.xB])
    # TD is the bubble point of saturated liquid reflux/distillate. The top
    # equilibrium tray temperature is the dew point of vapor y=xD.
    T_top_vapor=bubble_T(equilibrium_x(b.xD,b.P),b.P)
    alpha_avg=sqrt(relative_volatility(b.xD,b.P)*relative_volatility(b.xB,b.P))
    n_fenske=log(b.xD*(1-b.xB)/(b.xB*(1-b.xD)))/log(alpha_avg)
    Tref=(TD+TB)/2
    mu_eb,mu_sm=viscosity(Tref,"EB")*1000,viscosity(Tref,"SM")*1000
    # Logarithmic ideal mixture rule evaluated at feed composition; assumed.
    mu_mix=exp(b.zF*log(mu_eb)+(1-b.zF)*log(mu_sm))
    eta_o=.492*(mu_mix*alpha_avg)**(-.245)
    if not 0 < eta_o <= 1:
        raise ValueError("O'Connell correlation outside a physical range")
    actual=ceil((stages["N_fractional"]-1)/eta_o)
    murphree=murphree_steps(b,d,eta_o)
    sens=[]
    for factor in [1.1,1.2,1.5,2.,3.]:
        dd=operating_data(b,factor)
        ss=equilibrium_steps(b,dd)
        nn=packed_ntu(b,dd)
        sens.append(dict(R_over_Rmin=factor,R=dd["R"],V=dd["V"],
                         N_fractional=ss["N_fractional"],N_integer=ss["N_integer"],
                         feed_stage=ss["feed_stage"],NOG=nn["NOG_total"],
                         Z_at_assumed_HOG=nn["Z_total"]))
    grid=np.linspace(b.xB,b.xD,501)
    checks={
        "total_mass_residual_kmol_h": b.F-d["D"]-d["B"],
        "EB_balance_residual_kmol_h": b.F*b.zF-d["D"]*b.xD-d["B"]*b.xB,
        "q_line_intersection_residual":d["ar"]*b.zF+d["br"]-(d["ass"]*b.zF+d["bs"]),
        "minimum_operating_driving_force":min(equilibrium_y(x,b.P)-op_y(x,b,d) for x in grid),
        "maximum_bubble_pressure_residual_Pa":max(abs(x*psat(bubble_T(x,b.P),"EB")+(1-x)*psat(bubble_T(x,b.P),"SM")-b.P) for x in grid),
        "NTU_variable_change_error":ntu["integration_variable_check"],
        "minimum_count_lt_operating_count":minimum["N_fractional"]<stages["N_fractional"],
        "N_decreases_as_R_increases":all(sens[i]["N_fractional"]>sens[i+1]["N_fractional"] for i in range(len(sens)-1)),
    }
    limit=operating_data(b,1+1e-12)
    limit_grid=np.unique(np.r_[grid,b.zF])
    checks["minimum_reflux_other_pinches_absent"]=bool(min(equilibrium_y(x,b.P)-op_y(x,b,limit) for x in limit_grid)>=-1e-10)
    const=lambda x:alpha_avg*x/(1+(alpha_avg-1)*x)
    cs=equilibrium_steps(b,total_reflux=True,eq_override=const)
    checks["constant_alpha_integer_steps_equal_Fenske_ceiling"]=(cs["N_integer"]==ceil(n_fenske))
    assert abs(checks["total_mass_residual_kmol_h"])<1e-10
    assert abs(checks["EB_balance_residual_kmol_h"])<1e-10
    assert abs(checks["q_line_intersection_residual"])<1e-10
    assert checks["minimum_operating_driving_force"]>0
    assert checks["maximum_bubble_pressure_residual_Pa"]<1e-4
    assert checks["NTU_variable_change_error"]<1e-7
    assert all(checks[k] for k in ["minimum_count_lt_operating_count","N_decreases_as_R_increases","constant_alpha_integer_steps_equal_Fenske_ceiling","minimum_reflux_other_pinches_absent"])
    # Independent check: scaling F changes flows, not ideal N, Rmin or NTU.
    scaled=Basis(F=2*b.F,zF=b.zF,xD=b.xD,xB=b.xB,P=b.P,reflux_factor=b.reflux_factor,HOG=b.HOG)
    sd=operating_data(scaled)
    checks["flow_scale_invariance"]=(abs(sd["Rmin"]-d["Rmin"])<1e-12 and abs(sd["V"]-2*d["V"])<1e-9 and abs(sd["ar"]-d["ar"])<1e-12)
    assert checks["flow_scale_invariance"]
    return dict(basis=b.__dict__,flows=d,stages=stages,minimum=minimum,ntu=ntu,
                temperatures_K=dict(reflux_bubble=TD,feed=TF,bottom=TB,top_vapor_dew=T_top_vapor,viscosity_reference=Tref),
                alpha=dict(top=relative_volatility(b.xD,b.P),feed=relative_volatility(b.zF,b.P),bottom=relative_volatility(b.xB,b.P),average=alpha_avg),
                Nmin_Fenske_including_reboiler=n_fenske,
                viscosity_cP=dict(EB=mu_eb,SM=mu_sm,mixture=mu_mix),
                eta_overall=eta_o,actual_trays_OConnell=actual,murphree=murphree,
                sensitivity=sens,validation=checks)


def write_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def plot_results(b,r,out):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,
                         "axes.spines.top":False,"axes.spines.right":False,
                         "axes.grid":True,"grid.alpha":.18,"savefig.dpi":180})
    navy,teal,orange,gray="#18304b","#008d95","#db632f","#76879a"
    x=np.linspace(0,1,1001); y=np.array([equilibrium_y(xx,b.P) for xx in x]); d=r["flows"]
    T=np.array([bubble_T(xx,b.P)-273.15 for xx in x])
    def save(fig,name):
        fig.savefig(out/"images"/name,bbox_inches="tight",facecolor="white");plt.close(fig)
    def axes_xy(ax):
        ax.plot(x,y,color=navy,lw=2,label="Raoult equilibrium")
        ax.plot([0,1],[0,1],color=gray,lw=1,ls="--",label="$y=x$")
        ax.set(xlim=(0,1),ylim=(0,1),xlabel="Liquid EB mole fraction, x",ylabel="Vapor EB mole fraction, y")
        ax.set_aspect("equal",adjustable="box")
    fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout="constrained")
    axs[0].plot(x,T,color=navy,label="Bubble curve: T(x)")
    axs[0].plot(y,T,color=orange,label="Dew curve: T(y)")
    axs[0].set(xlabel="EB mole fraction (x or y)",ylabel="Temperature / deg C",title="T-x-y at 101.325 kPa")
    axs[0].legend();axes_xy(axs[1]);axs[1].set_title("Equilibrium separation is small")
    axs[1].legend(loc="upper left");save(fig,"01_vle.png")
    fig,axs=plt.subplots(1,2,figsize=(12,5.2),layout="constrained")
    for a in axs:axes_xy(a)
    p=np.array(r["minimum"]["points"]);axs[0].plot(p[:,0],p[:,1],color=orange,lw=1.0,label="Total-reflux steps")
    axs[0].set_title(f'Total reflux: Nmin = {r["minimum"]["N_fractional"]:.2f}\nIncludes one partial reboiler')
    axs[0].legend(loc="upper left",fontsize=9)
    rm=d["Rmin"];mr=rm/(rm+1);cr=b.xD/(rm+1)
    xs=np.linspace(b.xB,b.zF,150);xr=np.linspace(b.zF,b.xD,150)
    ms=(d["yp"]-b.xB)/(b.zF-b.xB)
    axs[1].plot(xr,mr*xr+cr,color=teal,label="Rectifying at Rmin")
    axs[1].plot(xs,ms*(xs-b.xB)+b.xB,color=orange,label="Stripping at Rmin")
    axs[1].axvline(b.zF,color=gray,ls=":",label="q = 1")
    axs[1].scatter([b.zF],[d["yp"]],color=orange,zorder=5)
    axs[1].annotate(f'Feed pinch\nRmin = {rm:.4f}',(b.zF,d["yp"]),xytext=(.13,.80),arrowprops=dict(arrowstyle="->",color=gray))
    axs[1].set_title("Minimum reflux: stage count diverges")
    axs[1].legend(loc="lower right",fontsize=8.5);save(fig,"02_limits.png")
    fig=plt.figure(figsize=(12,8),layout="constrained");gs=fig.add_gridspec(2,2,width_ratios=[1.45,1])
    ax=fig.add_subplot(gs[:,0]);axt=fig.add_subplot(gs[0,1]);axb=fig.add_subplot(gs[1,1])
    p=np.array(r["stages"]["points"])
    for a in [ax,axt,axb]:
        axes_xy(a);a.plot(xr,d["ar"]*xr+d["br"],color=teal,label="Rectifying")
        a.plot(xs,d["ass"]*xs+d["bs"],color=orange,label="Stripping")
        a.axvline(b.zF,color=gray,ls=":");a.plot(p[:,0],p[:,1],color="#a13741",lw=.85,label="Ideal stages")
    ax.set_title(f'McCabe-Thiele | R = {d["R"]:.4f}\n{r["stages"]["N_integer"]} stages including reboiler; feed stage {r["stages"]["feed_stage"]}')
    ax.legend(loc="upper left",fontsize=8.5)
    fr=r["stages"]["rows"][r["stages"]["feed_stage"]-1]
    ax.scatter([fr["x"]],[fr["y"]],s=35,color=orange,zorder=6)
    axt.set(xlim=(.88,1),ylim=(.88,1),title="Top-end detail")
    axb.set(xlim=(0,.10),ylim=(0,.10),title="Bottom-end detail")
    save(fig,"03_mccabe_thiele.png")
    fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout="constrained")
    for lo,hi,a,c,label in [(b.zF,b.xD,d["ar"],d["br"],"Rectifying"),(r["ntu"]["x_pack_bottom"],b.zF,d["ass"],d["bs"],"Stripping")]:
        xx=np.linspace(lo,hi,500); yy=a*xx+c;inv=1/np.array([equilibrium_y(v,b.P)-q for v,q in zip(xx,yy)])
        axs[0].plot(yy,inv,label=label,lw=2)
    axs[0].set(xlabel="Bulk vapor EB mole fraction, y",ylabel="$1/(y^*-y)$",title=f'Packing NTU = {r["ntu"]["NOG_total"]:.2f}\nPartial reboiler excluded')
    axs[0].legend()
    hh=np.linspace(.15,.85,200)
    axs[1].plot(hh,hh*r["ntu"]["NOG_total"],color=teal,lw=2)
    axs[1].scatter([b.HOG],[r["ntu"]["Z_total"]],color=orange,zorder=5)
    axs[1].annotate(f'Assumed HTU = {b.HOG:.2f} m\nPacking height = {r["ntu"]["Z_total"]:.2f} m',(b.HOG,r["ntu"]["Z_total"]),xytext=(.20,.88*max(hh)*r["ntu"]["NOG_total"]),arrowprops=dict(arrowstyle="->",color=gray))
    axs[1].set(xlabel="Assumed overall gas HTU / m",ylabel="Packing height / m",title="Height is conditional on HTU")
    save(fig,"04_ntu_htu.png")
    fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout="constrained")
    for rows,label,color in [(r["stages"]["rows"],"Ideal stages",teal),(r["murphree"]["rows"],"Murphree scenario",orange)]:
        nn=np.arange(1,len(rows)+1);xx=[row["x"] for row in rows]
        axs[0].plot(nn,xx,color=color,label=label,lw=2)
    axs[0].axhline(b.zF,color=gray,ls=":")
    axs[0].set(xlabel="Stage / tray number from top",ylabel="Liquid EB mole fraction",title=f'Overall efficiency estimate = {r["eta_overall"]:.3f}\nMurphree scenario assumes eta_M = eta_o')
    axs[0].legend()
    mus=np.linspace(.15,.35,200);es=.492*(mus*r["alpha"]["average"])**(-.245)
    axs[1].plot(mus,es,color=navy,lw=2)
    axs[1].scatter([r["viscosity_cP"]["mixture"]],[r["eta_overall"]],color=orange)
    axs[1].set(xlabel="Liquid mixture viscosity / cP",ylabel="Overall tray efficiency",title="O'Connell correlation, lecture 7 p.15")
    save(fig,"05_stage_efficiency.png")
    fig,axs=plt.subplots(1,2,figsize=(12,4.6),layout="constrained")
    ff=[v["R_over_Rmin"] for v in r["sensitivity"]]
    axs[0].plot(ff,[v["N_fractional"] for v in r["sensitivity"]],"o-",color=navy)
    axs[0].axhline(r["minimum"]["N_fractional"],color=gray,ls="--",label="Minimum stages")
    axs[0].set(xlabel="R / Rmin",ylabel="Theoretical stages, including reboiler",title="More reflux reduces stage count");axs[0].legend()
    axs[1].plot(ff,[v["V"] for v in r["sensitivity"]],"o-",color=orange)
    axs[1].set(xlabel="R / Rmin",ylabel="Internal vapor rate / (kmol/h)",title="More reflux increases internal traffic")
    save(fig,"06_reflux_sensitivity.png")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--feed",type=float,default=100.)
    ap.add_argument("--z-feed",type=float,default=.5)
    ap.add_argument("--x-top",type=float,default=.99)
    ap.add_argument("--x-bottom",type=float,default=.01)
    ap.add_argument("--reflux-factor",type=float,default=1.5)
    ap.add_argument("--htu",type=float,default=.5)
    ap.add_argument("--out",type=Path,default=Path(__file__).resolve().parent)
    args=ap.parse_args()
    b=Basis(F=args.feed,zF=args.z_feed,xD=args.x_top,xB=args.x_bottom,reflux_factor=args.reflux_factor,HOG=args.htu)
    out=args.out.resolve();(out/"images").mkdir(parents=True,exist_ok=True);(out/"data").mkdir(exist_ok=True)
    r=design(b)
    (out/"data"/"results.json").write_text(json.dumps(r,indent=2,ensure_ascii=False),encoding="utf-8")
    (out/"data"/"validation.json").write_text(json.dumps(r["validation"],indent=2),encoding="utf-8")
    write_csv(out/"data"/"ideal_stages.csv",r["stages"]["rows"])
    write_csv(out/"data"/"murphree_trays.csv",r["murphree"]["rows"])
    write_csv(out/"data"/"reflux_sensitivity.csv",r["sensitivity"])
    write_csv(out/"data"/"vle.csv",[dict(x=x,y=equilibrium_y(x,b.P),T_K=bubble_T(x,b.P),alpha=relative_volatility(x,b.P)) for x in np.linspace(0,1,1001)])
    plot_results(b,r,out)
    summary={k:r[k] for k in ["basis","flows","temperatures_K","alpha","Nmin_Fenske_including_reboiler","viscosity_cP","eta_overall","actual_trays_OConnell","ntu","validation"]}
    summary.update(N_ideal_total=r["stages"]["N_fractional"],N_ideal_total_integer=r["stages"]["N_integer"],feed_ideal=r["stages"]["feed_stage"],Nmin_variable_alpha=r["minimum"]["N_fractional"],feed_murphree_approx=r["murphree"]["feed_tray"])
    print(json.dumps(summary,indent=2))
    (out/"results_summary.md").write_text("# Recalculated results\n\n```json\n"+json.dumps(summary,indent=2)+"\n```\n\nREADME.md describes the original baseline; update its numbers if inputs change.\n",encoding="utf-8")


if __name__=="__main__":
    main()
