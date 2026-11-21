"""Task 21 -- manuscript figures 1-6 (Figure 6: RQ4 / H4, see fig6 docstring F6).

Plan: input/ideas/zero-shot-limits-context-complexity-plan.md, Task 21.  Verdicts shown follow knowledge/_meta/log.md:
H2a confirmed (amended A1), H2b rejected, H2c "holds for distinct cell types, not for donor lines", H1 untestable
(registered ceiling invalid: the oracle exceeds it), H3a confirmed (as registered), H3b inconclusive.

ASSUMPTIONS / RULES (stated before plotting; checked in code where possible)
  F1 Reads ONLY files in output/data-analysis/results/. No statistic is re-estimated: every point estimate and CI comes
     from the registered analysis outputs (rq2_h2a.json, rq2_params_tahoe.json, rq2_h2c.json, rq2_source_curve.json,
     rq1_ceiling.csv, rq1_h1.json, rq3_h3a.json, rq3_h3a_matchedn.json, rq3_h3b.json). Parquets are used only for
     per-run / per-subset display marks (per-seed lines, subset dots) and are cross-checked against the JSON points.
  F2 Exploratory quantities (include-zero curve, r_proj, matched-n D_other) are drawn with open marks / muted styling
     and labelled "exploratory" in the figure itself, never mixed into a registered mark.
  F3 Colours: the dataviz reference palette, slots 1-3 only (blue #2a78d6, orange #eb6834, aqua #1baf7a), validated
     all-pairs with the skill's validator (CVD dE 9.2, normal-vision 24.0; aqua < 3:1 contrast -> every aqua mark has a
     legend entry or direct label). References/ceilings use ink and greys. Every series also differs by marker shape.
  F4 Print: white surface, dark ink, vector PDF (TrueType fonts embedded), 6.5 in wide, fonts >= 8 pt.
  F5 The only randomness is the horizontal jitter of subset dots in Fig 3b (seed 20261121).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
OUT = HERE.parent / "article" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261121
W = 6.5  # inches, single-column arXiv text width

# ---- palette (dataviz reference instance, light/print) ----------------------------------------------------------
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
BLUE_TINT, GRAY_TINT = "#cde2fb", "#f0efec"
ORANGE_TINT = "#fbe1d6"   # S2 at ~20 % over white (wash, not a new hue)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.edgecolor": AXIS, "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "axes.labelcolor": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "grid.color": GRID, "grid.linewidth": 0.5, "grid.linestyle": "-",
    "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round", "legend.frameon": False,
    "axes.titlelocation": "left", "axes.titleweight": "normal", "mathtext.fontset": "dejavusans",
})
LW = 1.5          # ~2 px line
MS = 5.5          # marker size in pt (>= 8 px)
PANELS = ["H1", "OR_HCT116X", "KOLF", "RN_K562", "RN_HepG2", "RN_Jurkat", "RN_RPE1"]
PLAB = {"H1": "H1", "OR_HCT116X": "HCT116 (OR)", "KOLF": "KOLF2.1J", "RN_K562": "K562 (RN)",
        "RN_HepG2": "HepG2 (RN)", "RN_Jurkat": "Jurkat (RN)", "RN_RPE1": "RPE1 (RN)"}


def jload(name):
    return json.loads((RES / name).read_text())


def tag(ax, letter, x=-0.02, y=1.02):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=220)
    plt.close(fig)
    print("wrote", OUT / f"{name}.pdf")


def ci_h(ax, y, x, lo, hi, color, marker="o", filled=True, ms=MS, lw=1.0, z=3, label=None):
    ax.plot([lo, hi], [y, y], color=color, lw=lw, zorder=z, solid_capstyle="butt")
    ax.plot([x], [y], marker=marker, ms=ms, color=color, mfc=color if filled else "white", mec=color,
            mew=1.1, ls="none", zorder=z + 1, label=label)


# ======================================================================================================== Fig 1 =====
def fig1():
    fig = plt.figure(figsize=(W, 3.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 56)
    ax.axis("off")

    def box(x, y, w, h, fc, ec=AXIS, lw=0.7, r=1.2):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw))

    def arrow(p, q, color=INK2, lw=1.0):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=9, color=color, lw=lw,
                                     shrinkA=0, shrinkB=0))

    def pill(x, y, text):
        ax.text(x, y, text, fontsize=8, fontweight="bold", color="white", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.6", fc=INK, ec=INK))

    # --- top: decomposition ---------------------------------------------------------------------------------
    ax.text(1.5, 54.3, "Response decomposition (pseudobulk; context c, target t)", fontsize=8.5, va="center")
    yb, hb, wb = 39.5, 6.5, 6.5
    xs = [10.0, 19.5, 29.0, 38.5, 48.0]
    terms = [(r"$\mu$", BLUE_TINT, "shared\nresponse"), (r"$a_t$", BLUE_TINT, "target\neffect"),
             (r"$b_c$", ORANGE_TINT, "context\neffect"), (r"$i_{c,t}$", ORANGE_TINT, "context ×\ntarget"),
             (r"$\varepsilon_{c,t}$", GRAY_TINT, "noise")]
    ax.text(4.0, yb + hb / 2, r"$\delta_{c,t}$", fontsize=11, ha="center", va="center")
    ax.text(8.0, yb + hb / 2, "=", fontsize=11, ha="center", va="center")
    for k, (x, (sym, fc, lab)) in enumerate(zip(xs, terms)):
        box(x, yb, wb, hb, fc)
        ax.text(x + wb / 2, yb + hb / 2, sym, fontsize=11, ha="center", va="center")
        ax.text(x + wb / 2, yb - 0.9, lab, fontsize=8, ha="center", va="top", color=INK2, linespacing=1.05)
        if k:
            ax.text(x - 1.5, yb + hb / 2, "+", fontsize=11, ha="center", va="center")
    lx = 62.0
    for j, (fc, t) in enumerate([(BLUE_TINT, "transferable from source contexts"),
                                 (ORANGE_TINT, "context-specific: predicted from $x_{c^*}$"),
                                 (GRAY_TINT, "noise: removed by split halves")]):
        box(lx, 44.0 - j * 3.9, 2.6, 2.4, fc, r=0.5)
        ax.text(lx + 3.8, 45.2 - j * 3.9, t, fontsize=8, va="center")
    xb0, xb1, yr = xs[0], xs[2] + wb, 48.0
    ax.plot([xb0, xb0, xb1, xb1], [yr - 0.9, yr, yr, yr - 0.9], color=INK, lw=0.8)
    ax.plot([(xb0 + xb1) / 2] * 2, [yr, yr + 0.9], color=INK, lw=0.8)
    pill(12.5, 50.6, "RQ1")
    ax.text(15.6, 50.6, r"ceilings: $r_{\max}=\sqrt{(\mathrm{Var}\,a+\mathrm{Var}\,b)\,/\,\mathrm{Var}_{\mathrm{rep}}}$"
            r", oracle $\mu+a_t+b_c$, measurement ceiling", fontsize=8, va="center")

    ax.plot([1.5, 98.5], [31.5, 31.5], color=GRID, lw=0.8)
    ax.text(1.5, 29.3, "Zero-shot setting", fontsize=8.5, va="center")

    # --- bottom: zero-shot pipeline ---------------------------------------------------------------------------
    for j in (2, 1):
        box(2.0 + j * 0.9, 9.0 + j * 1.5, 24.0, 11.0, "white", ec=AXIS)
    box(2.0, 9.0, 24.0, 11.0, "white", ec=INK2)
    ax.text(3.3, 17.6, "C source contexts (c = 1 … C)", fontsize=8, va="center")
    ax.text(3.3, 14.4, "controls → descriptor $x_c$", fontsize=8, va="center", color=INK2)
    ax.text(3.3, 11.3, r"responses $\delta_{c,t}$ (many t)", fontsize=8, va="center", color=INK2)
    ax.text(2.0, 4.6, "C = independent samples for the map\n(not cells, not targets)", fontsize=8, va="center",
            color=INK2, linespacing=1.1)

    box(30.0, 10.0, 30.0, 11.0, BLUE_TINT, ec=S1, lw=0.9)
    ax.text(45.0, 18.3, "learned context map", fontsize=8, ha="center", va="center")
    ax.text(45.0, 15.0, r"$x\ \mapsto\ (\hat b_c,\ \hat\imath_{c,t})$", fontsize=9, ha="center", va="center")
    ax.text(45.0, 11.9, r"identifiable rank $\leq C-1$ (Prop. 1)", fontsize=8, ha="center", va="center",
            color=INK2)
    arrow((26.2, 15.0), (29.8, 15.0))
    pill(32.5, 24.3, "RQ2")
    ax.text(35.6, 24.3, "context count: gain vs C, break-even $C^*$", fontsize=8, va="center")

    box(70.0, 18.0, 27.0, 7.0, "white", ec=INK2)
    ax.text(83.5, 23.0, "unseen context $c^*$", fontsize=8, ha="center", va="center")
    ax.text(83.5, 19.9, "control cells only → $x_{c^*}$", fontsize=8, ha="center", va="center", color=INK2)
    arrow((69.8, 21.0), (60.2, 18.0))

    box(64.0, 6.5, 34.0, 8.0, ORANGE_TINT, ec=S2, lw=0.9)
    ax.text(81.0, 12.4, "prediction for every target t", fontsize=8, ha="center", va="center")
    ax.text(81.0, 9.0, r"$\hat\delta_{c^*,t}=\hat\mu+\hat a_t+\hat b(x_{c^*})+\hat\imath_t(x_{c^*})$", fontsize=9,
            ha="center", va="center")
    arrow((60.2, 13.0), (63.8, 10.5))
    pill(66.5, 3.2, "RQ3")
    ax.text(69.6, 3.2, "residual subspace: does the error\nlie on $c^*$'s own control PCs?", fontsize=8,
            va="center", linespacing=1.1)
    save(fig, "fig1")

# ======================================================================================================== Fig 2 =====
def fig2():
    h2a = jload("rq2_h2a.json")
    par = jload("rq2_params_tahoe.json")["split_half"]["estimates"]
    conf = h2a["results"]["confirmatory"]
    expl = h2a["results"]["exploratory_include_zero"]
    prim = conf["summary"]["ctx_pathway_gain"]
    regq = conf["summary"]["mse_gain"]
    izq = expl["summary"]["ctx_pathway_gain"]
    Cs = np.array([r["C"] for r in prim["curve"]], float)
    sc = 1e4

    # per-seed lines from the parquet (frozen def.: unweighted mean over held-out contexts, then over folds)
    df = pd.read_parquet(RES / "rq2_tahoe_curve.parquet")
    ps = df.groupby(["C_requested", "seed", "fold"]).ctx_pathway_gain.mean().groupby(["C_requested", "seed"]).mean()
    seeds = sorted(df.seed.unique())
    for s in seeds:  # cross-check against the registered per-seed curves
        js = np.array(conf["per_seed"]["ctx_pathway_gain"]["seeds"][str(s)]["curve"])
        assert np.allclose(ps.xs(s, level="seed").values, js, atol=1e-12), f"per-seed mismatch seed {s}"
    mean_pq = ps.groupby("C_requested").mean().values
    assert np.allclose(mean_pq, [r["gain"] for r in prim["curve"]], atol=1e-12), "run-mean mismatch"

    fig = plt.figure(figsize=(W, 3.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.75, 1.0], left=0.105, right=0.985, bottom=0.34, top=0.91,
                          wspace=0.28)
    ax = fig.add_subplot(gs[0])
    bx = fig.add_subplot(gs[1], sharey=ax)

    def curve(a, q, color, z=3, band=True, label=None, marker="o", filled=True):
        g = np.array([r["gain"] for r in q["curve"]]) * sc
        lo = np.array([r["ci95"][0] for r in q["curve"]]) * sc
        hi = np.array([r["ci95"][1] for r in q["curve"]]) * sc
        if band:
            a.fill_between(Cs, lo, hi, color=color, alpha=0.13, lw=0, zorder=z - 1)
        a.plot(Cs, g, color=color, lw=LW, zorder=z, label=label, marker=marker, ms=MS - 1,
               mfc=color if filled else "white", mec="white" if filled else color, mew=1.0 if filled else 1.1)

    for a in (ax, bx):
        a.set_xscale("log")
        a.set_xticks(Cs)
        a.set_xticklabels([str(int(c)) for c in Cs])
        a.minorticks_off()
        a.axhline(0, color=INK2, lw=0.7, zorder=1)
        a.grid(axis="y")
        a.set_xlim(5.4, 47)
        a.set_xlabel("training contexts C (log scale)")

    for s in seeds:
        ax.plot(Cs, ps.xs(s, level="seed").values * sc, color=S1, lw=0.7, alpha=0.45, zorder=2)
    curve(ax, prim, S1, label="run mean, 95 % cluster-bootstrap CI")
    ax.set_ylabel("context-pathway gain\n(residual MSE, ×10$^{-4}$)")
    y0, y1 = -8.2, 3.6
    ax.set_ylim(y0, y1)
    fy = lambda v: (v - y0) / (y1 - y0)
    # C*_emp with CI (upper bound censored > 40)
    ce, (clo, _) = prim["cstar_emp"], prim["cstar_emp_ci95"]
    ax.plot([clo, 46.5], [0, 0], color=INK, lw=2.2, zorder=2.5, solid_capstyle="butt")
    ax.annotate("", xy=(47.0, 0), xytext=(45.0, 0), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0,
                                                                    mutation_scale=7), zorder=6)
    ax.plot([ce], [0], marker="D", ms=MS, color=INK, mec="white", mew=1.0, zorder=7)
    ax.annotate(f"$C^*_{{emp}}$ = {ce:.1f} [{clo:.1f}, > 40]", xy=(ce, -0.15), xytext=(45.5, -2.6), fontsize=8,
                ha="right", va="center", arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6, shrinkA=2))
    # off-axis theory values
    cmin, cmin_ci = par["C_min"]["estimate"], par["C_min"]["ci95"]
    cth, cth_ci = par["C_theory"]["estimate"], par["C_theory"]["ci95"]
    ax.text(0.02, fy(3.0), "Prop. 2 (rejected): break-evens beyond the axis", transform=ax.transAxes, fontsize=8,
            ha="left", va="center", color=INK2, style="italic")
    for yy, txt in ((2.05, f"$C^*_{{min}}$ = {cmin:.1f} [{cmin_ci[0]:.1f}, {cmin_ci[1]:.1f}]"),
                    (1.05, f"$C^*_{{theory}}$ = {cth:.0f} [{cth_ci[0]:.0f}, {cth_ci[1]:.0f}]")):
        ax.text(0.02, fy(yy), txt, transform=ax.transAxes, fontsize=8, ha="left", va="center")
        ax.annotate("", xy=(1.0, fy(yy)), xycoords="axes fraction", xytext=(0.62, fy(yy)),
                    textcoords="axes fraction", arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.8,
                                                                mutation_scale=8))
    ax.set_title("Tahoe: context pathway vs constant-z ablation (primary, A1)", fontsize=8.5)
    tag(ax, "a", x=-0.13)
    h_ps = Line2D([], [], color=S1, lw=0.7, alpha=0.6, label="per-seed mean (3 seeds)")
    h_m = Line2D([], [], color=S1, lw=LW, marker="o", ms=MS - 1, mec="white", label="run mean ± 95 % CI")
    h_c = Line2D([], [], color=INK, marker="D", ms=MS, mec="white", lw=2.0, label="$C^*_{emp}$ with 95 % CI")
    ax.legend(handles=[h_m, h_ps, h_c], loc="upper left", bbox_to_anchor=(-0.02, -0.24), ncol=2, fontsize=8,
              handlelength=1.6, columnspacing=1.0, labelspacing=0.3, borderaxespad=0)

    curve(bx, prim, S1, z=5, band=False, label=f"primary (A1), $C^*$ {prim['cstar_emp']:.1f}")
    curve(bx, regq, S2, z=4, label=f"registered (confounded), $C^*$ {regq['cstar_emp']:.1f}", marker="s")
    curve(bx, izq, S3, z=3, label=f"include-zero (exploratory), $C^*$ {izq['cstar_emp']:.1f}",
          marker="^", filled=False)
    plt.setp(bx.get_yticklabels(), visible=False)
    bx.set_title("same runs, other quantities", fontsize=8.5)
    tag(bx, "b", x=-0.04)
    bx.legend(loc="upper left", bbox_to_anchor=(-0.12, -0.24), ncol=1, fontsize=8, handlelength=1.6,
              borderaxespad=0, labelspacing=0.3)
    save(fig, "fig2")


# ======================================================================================================== Fig 3 =====
def fig3():
    h2c = jload("rq2_h2c.json")
    par = jload("rq2_params_tahoe.json")["split_half"]["estimates"]
    h2a = jload("rq2_h2a.json")["results"]["confirmatory"]["summary"]["ctx_pathway_gain"]
    sc = jload("rq2_source_curve.json")
    gene = jload("rq2_params_gene.json")
    c1, c2 = h2c["C_public_count1"], h2c["C_public_count2"]
    assert (c1, c2) == (gene["count2"]["C_public_count1"], gene["count2"]["C_public_count2"])
    lb = h2c["S1_S2_S3"]["line_bootstrap"]["C_min"]["ci95"]

    fig = plt.figure(figsize=(W, 5.6))
    top = fig.add_gridspec(1, 1, left=0.335, right=0.975, top=0.925, bottom=0.675)
    ax = fig.add_subplot(top[0])
    rows = [  # label, point, lo, hi, color, marker, filled, censored
        ("Genetic $C^*_{min}$ (registered)", h2c["C_min"], *h2c["C_min_ci95"], S1, "o", True, False),
        ("Genetic $C^*_{min}$, line bootstrap (S3)", h2c["C_min"], *lb, S1, "o", False, False),
        ("Tahoe $C^*_{emp}$ (measured)", h2a["cstar_emp"], h2a["cstar_emp_ci95"][0], 40.0, S2, "D", True,
         True),
        ("Tahoe $C^*_{min}$ (same estimator)", par["C_min"]["estimate"], *par["C_min"]["ci95"], S2, "o", False, False),
        ("Tahoe $C^*_{theory}$ (Prop. 2; rejected)", par["C_theory"]["estimate"], *par["C_theory"]["ci95"], S2, "s",
         False, False),
    ]
    ys = np.arange(len(rows))[::-1].astype(float)
    ys[2:] -= 0.35
    for c, lab in ((c1, f"Count 1 = {c1}\ncell types"), (c2, f"Count 2 = {c2}\nincl. donor lines")):
        ax.axvline(c, color=INK2, lw=0.8, zorder=1)
        ax.text(c * 1.04, ys[0] + 1.75, lab, fontsize=8, ha="left", va="top", color=INK2, linespacing=1.0)
    for y, (lab, x, lo, hi, col, mk, fill, cens) in zip(ys, rows):
        ci_h(ax, y, x, lo, hi, col, marker=mk, filled=fill, lw=1.6 if fill else 1.0)
        if cens:
            ax.annotate("", xy=(47, y), xytext=(39.0, y), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0,
                                                                          mutation_scale=7))
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_xscale("log")
    ax.set_xticks([5, 10, 20, 50, 100, 200, 300])
    ax.set_xticklabels(["5", "10", "20", "50", "100", "200", "300"])
    ax.minorticks_off()
    ax.set_xlim(5, 330)
    ax.set_ylim(ys[-1] - 0.6, ys[0] + 1.9)
    ax.grid(axis="x")
    ax.set_xlabel("number of contexts (log scale)")
    ax.spines["left"].set_visible(False)
    # calibration gap bracket
    ratio = par["C_min"]["estimate"] / h2a["cstar_emp"]
    yg = (ys[2] + ys[3]) / 2
    ax.annotate("", xy=(h2a["cstar_emp"], yg), xytext=(par["C_min"]["estimate"], yg),
                arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.7, mutation_scale=7))
    ax.text(np.sqrt(h2a["cstar_emp"] * par["C_min"]["estimate"]), yg + 0.08, f"×{ratio:.1f}", fontsize=8,
            ha="center", va="bottom", color=INK2)
    ax.text(320, (ys[0] + ys[1]) / 2, f"H2c: {h2c['verdict']}", fontsize=8, ha="right", va="center",
            style="italic", color=INK)
    ax.set_title("Break-even number of contexts vs public genetic data", fontsize=8.5)
    tag(ax, "a", x=-0.47, y=1.02)

    # ---- panel b: source-count curves -------------------------------------------------------------------------------
    df = pd.read_parquet(RES / "rq2_source_curve.parquet")
    df = df[df.arm == "own"]
    rng = np.random.default_rng(SEED)
    bot = fig.add_gridspec(2, 4, left=0.085, right=0.985, top=0.50, bottom=0.075, wspace=0.50, hspace=0.62)
    axes = [fig.add_subplot(bot[i // 4, i % 4]) for i in range(8)]
    for i, p in enumerate(PANELS):
        a = axes[i]
        d = df[df.panel == p]
        K = int(d.K.iloc[0])
        cur = sc["panels"][p]["own"]["curve"]
        ks = np.array(sorted(int(k) for k in cur))
        mu = np.array([cur[str(k)]["pds"] for k in ks])
        assert np.allclose(d.groupby("k").pds.mean().reindex(ks).values, mu, atol=1e-10), p
        jit = np.exp(rng.uniform(-0.07, 0.07, len(d)))
        a.scatter(d.k.values * jit, d.pds.values, s=4, color=MUTED, alpha=0.45, lw=0, zorder=2)
        a.plot(ks, mu, color=S1, lw=LW, marker="o", ms=MS - 1.5, mec="white", mew=0.8, zorder=4)
        a.set_xscale("log")
        a.set_xticks([1, 2, 4, K])
        a.set_xticklabels(["1", "2", "4", str(K)])
        a.minorticks_off()
        a.set_xlim(0.85, K * 1.15)
        a.grid(axis="y")
        a.set_title(PLAB[p], fontsize=8)
        a.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
        a.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.2f"))
        if i % 4 == 0:
            a.set_ylabel("PDS (offline)")
        if i >= 4:
            a.set_xlabel("source families k")
    tag(axes[0], "b", x=-0.32, y=1.22)
    # LOFO forest
    a = axes[7]
    yy = np.arange(len(PANELS))[::-1]
    for y, p in zip(yy, PANELS):
        lf = sc["panels"][p]["own"]["lofo"]
        ci_h(a, y, lf["mean_dPDS"] * 100, lf["ci95"][0] * 100, lf["ci95"][1] * 100, S1, ms=MS - 1)
    a.axvline(0, color=INK2, lw=0.7, zorder=1)
    a.set_yticks(yy)
    a.set_yticklabels([PLAB[p].replace(" (RN)", "").replace(" (OR)", "").replace("2.1J", "") for p in PANELS])
    a.tick_params(axis="y", length=0)
    a.spines["left"].set_visible(False)
    a.grid(axis="x")
    a.set_xlabel("LOFO ΔPDS (×10$^{-2}$)")
    a.set_title("LOFO marginal value", fontsize=8)
    h = [Line2D([], [], color=MUTED, marker="o", ls="none", ms=3, alpha=0.6, label="family subset"),
         Line2D([], [], color=S1, lw=LW, marker="o", ms=MS - 1.5, mec="white", label="mean over subsets")]
    fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.075, 0.585), ncol=2, fontsize=8, handlelength=1.6)
    save(fig, "fig3")


# ======================================================================================================== Fig 4 =====
def fig4():
    d = pd.read_csv(RES / "rq1_ceiling.csv")
    h1 = jload("rq1_h1.json")["H1"]
    comp = jload("rq1_components.json")

    def get(p, pred, q):
        r = d[(d.panel == p) & (d.predictor == pred) & (d.quantity == q)]
        if r.empty:
            return np.nan, np.nan, np.nan
        r = r.iloc[0]
        return r.estimate, r.ci_lo, r.ci_hi

    assert abs(get("H1", "ceil", "r_max")[0] - h1["r_max"]["est"]) < 1e-12
    assert abs(get("H1", "m3", "r")[0] - h1["r_ours_m3"]["est"]) < 1e-12
    assert "H1" in comp["panels"]
    oracle_pd = h1["validity"]["per_predictor"]["oracle"]["paired_diff_ci"]

    fig = plt.figure(figsize=(W, 3.75))
    ax = fig.add_axes([0.13, 0.115, 0.84, 0.70])
    ys = np.arange(len(PANELS))[::-1].astype(float) * 1.0
    ys[0] += 0.35   # room for the H1 annotation
    up, dn = 0.17, -0.17
    for y, p in zip(ys, PANELS):
        ax.axhspan(y - 0.42, y + 0.42, color=GRAY_TINT if p == "H1" else "white", lw=0, zorder=0)
        cm = get(p, "ceil", "ceiling_meas")
        ax.plot([cm[0]] * 2, [y - 0.32, y + 0.32], color=INK, lw=1.8, zorder=3, solid_capstyle="butt")
        rmx = get(p, "ceil", "r_max")
        if np.isfinite(rmx[0]):
            ci_h(ax, y + up, rmx[0], rmx[1], rmx[2], INK2, marker="D", filled=False, ms=MS)
        else:
            ax.text(0.40, y + up, "$r_{max}$ and oracle not estimable", fontsize=8, va="center", color=INK2,
                    style="italic")
        rp = get(p, "ceil", "r_proj")
        if np.isfinite(rp[0]):
            ci_h(ax, y + up, rp[0], rp[1], rp[2], MUTED, marker="s", filled=False, ms=MS - 0.5, lw=0.8)
        orc = get(p, "oracle", "r")
        if np.isfinite(orc[0]):
            ci_h(ax, y + up, orc[0], orc[1], orc[2], INK, marker="o", filled=True, ms=MS - 2, lw=0.0, z=6)
        m3 = get(p, "m3", "r")
        ci_h(ax, y + dn, m3[0], m3[1], m3[2], S1, marker="o", filled=True, lw=1.4)
        xm = get(p, "xmean", "r")
        ci_h(ax, y + dn, xm[0], xm[1], xm[2], S2, marker="^", filled=True, lw=1.4)
    # H1 annotation
    yH = ys[0]
    rmx = get("H1", "ceil", "r_max")[0]
    ax.annotate("registered $r_{max}$ (H1): invalid, the oracle exceeds it\n"
                f"(paired diff. [{oracle_pd[0]:+.3f}, {oracle_pd[1]:+.3f}])  →  H1 untestable",
                xy=(rmx, yH + up + 0.06), xytext=(0.03, yH + 0.80), fontsize=8, ha="left", va="center",
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.7, mutation_scale=7, shrinkB=3))
    ax.set_yticks(ys)
    ax.set_yticklabels([PLAB[p] for p in PANELS])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(ys[-1] - 0.55, ys[0] + 1.2)
    ax.grid(axis="x")
    ax.set_xlabel("noise-corrected Pearson r with the destination response (target bootstrap 95 % CI)")
    h = [Line2D([], [], color=INK, lw=1.8, marker="|", ms=9, ls="none", mew=1.8, label="measurement ceiling"),
         Line2D([], [], color=INK2, marker="D", mfc="white", ms=MS, lw=1.0, label="$r_{max}$ (registered estimator)"),
         Line2D([], [], color=INK, marker="o", ms=MS - 2, ls="none", label="oracle $\\mu+a_t+b_c$ (diagnostic)"),
         Line2D([], [], color=MUTED, marker="s", mfc="white", ms=MS - 0.5, lw=0.8, label="$r_{proj}$ (exploratory)"),
         Line2D([], [], color=S1, marker="o", ms=MS, lw=1.4, label="pipeline M1+M2+M3"),
         Line2D([], [], color=S2, marker="^", ms=MS, lw=1.4, label="cross-line mean")]
    fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.02, 0.995), ncol=3, fontsize=8, handlelength=1.8,
               columnspacing=1.2, labelspacing=0.35)
    ax.text(1.0, 1.01, "upper lane: ceilings · lower lane: zero-shot methods", transform=ax.transAxes, fontsize=8,
            ha="right", va="bottom", color=INK2)
    save(fig, "fig4")


# ======================================================================================================== Fig 5 =====
def fig5():
    a3 = jload("rq3_h3a.json")
    mn = jload("rq3_h3a_matchedn.json")
    b3 = jload("rq3_h3b.json")
    assert mn["label"].startswith("exploratory")
    for p in PANELS:  # matched-n file reproduces the registered numbers
        assert abs(mn["table"][p]["k8"]["reg"]["D_other"] - a3["panels"][p]["k8"]["D_other"]["point"]) < 1e-9, p

    fig = plt.figure(figsize=(W, 5.0))
    gtop = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], left=0.14, right=0.985, top=0.88, bottom=0.53,
                            wspace=0.10)
    ax = fig.add_subplot(gtop[0])
    bx = fig.add_subplot(gtop[1], sharey=ax)
    ys = np.arange(len(PANELS))[::-1].astype(float)
    sub = [("F_dest", "destination top-8 PCs", S1, "o", True),
           ("F_other_mean", "other lines' PCs (mean of 6)", S2, "s", True),
           ("F_perm_mean", "permuted-control PCs", S3, "^", False),
           ("F_haar", "random (Haar) 8-subspace", MUTED, "D", False)]
    off = {"F_dest": 0.21, "F_other_mean": 0.07, "F_perm_mean": -0.07, "F_haar": -0.21}
    for y, p in zip(ys, PANELS):
        k8 = a3["panels"][p]["k8"]
        if p == "H1":
            ax.axhspan(y - 0.45, y + 0.45, color=GRAY_TINT, lw=0, zorder=0)
            bx.axhspan(y - 0.45, y + 0.45, color=GRAY_TINT, lw=0, zorder=0)
        for q, lab, col, mk, fill in sub:
            v = k8[q]
            ci_h(ax, y + off[q], v["point"], v["target_ci"][0], v["target_ci"][1], col, marker=mk, filled=fill,
                 ms=MS - 0.5, lw=1.0, label=lab if p == "H1" else None)
        reg = k8["D_other"]
        m = mn["table"][p]["k8"]["mn"]
        ci_h(bx, y + 0.12, reg["point"], *reg["target_ci"], S1, marker="o", filled=True, lw=1.3,
             label="registered" if p == "H1" else None)
        ci_h(bx, y - 0.12, m["D_other"], *m["target_ci"], S2, marker="o", filled=False, lw=1.0,
             label="matched n (exploratory)" if p == "H1" else None)
    ax.set_xscale("log")
    ax.set_xlim(8e-4, 0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([PLAB[p] for p in PANELS])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x")
    ax.set_xlabel("residual energy fraction F(U), k = 8 (log scale)")
    ax.set_title("H3a: where the transfer residual's energy lies", fontsize=8.5)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.2, 1.07), ncol=2, fontsize=8, handlelength=1.4,
              columnspacing=1.0, labelspacing=0.3)
    tag(ax, "a", x=-0.20, y=1.25)
    bx.axvline(0, color=INK2, lw=0.7, zorder=1)
    plt.setp(bx.get_yticklabels(), visible=False)
    bx.tick_params(axis="y", length=0)
    bx.spines["left"].set_visible(False)
    bx.grid(axis="x")
    bx.set_xlabel("$D_{other}$ = F(U$_{dest}$) − mean F(U$_j$)")
    bx.set_title("destination vs other lines", fontsize=8.5)
    bx.legend(loc="lower left", bbox_to_anchor=(0.0, 1.07), ncol=1, fontsize=8, handlelength=1.4,
              labelspacing=0.3)
    tag(bx, "b", x=-0.04, y=1.25)

    # ---- panel c: H3b --------------------------------------------------------------------------------------------
    cx = fig.add_axes([0.25, 0.09, 0.735, 0.31])
    alts = [("none", "no correction"), ("M3_other", "other-line PC deflation"), ("M3_perm", "permuted-PC deflation"),
            ("predPC", "prediction-PC deflation"), ("sys_best", "Systema-style (best η)")]
    real = b3["real_scorer_H1_dPDS"]
    yc = np.arange(len(alts))[::-1].astype(float)
    for y, (k, lab) in zip(yc, alts):
        v = b3["dPDS_7panel"][k]
        pp = v["per_panel"]
        others = [pp[p] for p in PANELS if p != "H1"]
        cx.scatter(np.array(others) * 100, np.full(len(others), y - 0.2), s=12, color=MUTED, alpha=0.7, lw=0,
                   zorder=2, label="per panel (6 non-H1)" if k == "none" else None)
        cx.plot([pp["H1"] * 100], [y - 0.2], marker="o", ms=MS - 1, mfc="white", mec=INK, mew=1.0, ls="none",
                zorder=4, label="H1, offline" if k == "none" else None)
        ci_h(cx, y + 0.12, v["mean"] * 100, v["ci95"][0] * 100, v["ci95"][1] * 100, S1, marker="o", lw=1.6,
             label="7-panel mean, 95 % CI" if k == "none" else None)
        if k in real:
            cx.plot([real[k]["dPDS_real"] * 100], [y - 0.36], marker="D", ms=MS, color=S2, mec="white", mew=0.8,
                    ls="none", zorder=5, label="H1, real scorer" if k == "none" else None)
    cx.axvline(0, color=INK2, lw=0.7, zorder=1)
    cx.set_yticks(yc)
    cx.set_yticklabels(["M3 vs " + a[1] for a in alts])
    cx.tick_params(axis="y", length=0)
    cx.spines["left"].set_visible(False)
    cx.grid(axis="x")
    cx.set_xlabel("ΔPDS = PDS(M3, destination PCs) − PDS(alternative)   (×10$^{-2}$)")
    cx.set_ylim(yc[-1] - 0.6, yc[0] + 0.55)
    cx.set_title(f"H3b: correction gain against five alternatives — verdict {b3['verdict']}", fontsize=8.5)
    cx.set_xlim(right=4.4)
    cx.legend(loc="upper right", fontsize=8, handlelength=1.2, labelspacing=0.3, borderaxespad=0.2)
    tag(cx, "c", x=-0.31, y=1.03)
    save(fig, "fig5")


# ======================================================================================================== Fig 6 =====
def fig6():
    """RQ4 / H4: same-line data value vs modelling gains (H1 panel, real scorer unless stated).

    F6 Data: rq4_h4.json (H1 real-scorer composites; gains = composite - composite(w = 0), recomputed and asserted
       against the stored gains), rq2_arch/rq2_arch_seeds_summary.json (best-architecture gain on the K562/HepG2
       screen: additive, two-line seed mean with t-CI over 3 seeds; largest single-line mean = K562 sensitivity),
       rq2_source_curve.json (H1 LOFO per-context value, offline sum-side composite, bootstrap CI). Linear axis so
       that negative CI ends stay visible (a log axis would hide them). The H1 runs are single runs; the emission-seed
       range (rq4_h4.json reference_emission_seed_spread) is drawn as a band, not as a CI.
    """
    h4 = jload("rq4_h4.json")
    arch = jload("rq2_arch/rq2_arch_seeds_summary.json")["per_arch"]
    lofo = jload("rq2_source_curve.json")["panels"]["H1"]["own"]["lofo"]
    c, g = h4["composites"], h4["gains"]
    ws = [0.0, 0.25, 0.5, 0.75, 1.0]
    gw = np.array([c[f"w{w:g}"] - c["w0"] for w in ws])
    for w, v in zip(ws[1:], gw[1:]):
        assert abs(v - g[f"H1GW_w{w:g}"]) < 1e-12, w
    gM2, gM3 = c["w0"] - c["noM2"], c["w0"] - c["noM3"]
    assert abs(gM2 - g["M2"]) < 1e-12 and abs(gM3 - g["M3"]) < 1e-12
    # best learned architecture = largest two-line seed mean over the six re-seeded architectures
    best = max(arch, key=lambda a: arch[a]["two_line_mean"]["mean"])
    tl = arch[best]["two_line_mean"]
    sens_line = max(("K562", "HepG2"), key=lambda L: arch[best][L]["mean"])
    sens = arch[best][sens_line]["mean"]
    assert best == "additive" and tl["n"] == 3 and sens_line == "K562"
    spread = h4["reference_emission_seed_spread"]["range"]
    largest = max(gM2, gM3, tl["mean"])
    ratio = gw[-1] / largest
    assert h4["verdict"] == "confirmed" and gw[-1] > largest
    pcv, pcv_ci = lofo["mean_dcomposite_sum"], lofo["ci95_composite_sum"]

    fig = plt.figure(figsize=(W, 3.35))
    # ---- panel a: gain vs H1GW weight ----------------------------------------------------------------------------
    ax = fig.add_axes([0.075, 0.20, 0.255, 0.61])
    ax.fill_between(ws, gw - spread / 2, gw + spread / 2, color=ORANGE_TINT, lw=0, zorder=1)
    ax.plot(ws, gw, color=S2, lw=LW, marker="o", ms=MS, mec="white", mew=0.8, zorder=4)
    ax.axhline(0, color=INK2, lw=0.7, zorder=2)
    ax.axhline(tl["mean"], color=S1, lw=1.0, ls=(0, (4, 2)), zorder=2)
    ax.text(1.02, tl["mean"] + 0.0012, f"largest modelling gain\n(arch., +{tl['mean']:.4f})", fontsize=8, color=INK2,
            ha="right", va="bottom", linespacing=1.0)
    ax.annotate(f"w = 1 (registered)\n+{gw[-1]:.4f}", xy=(1.0, gw[-1]), xytext=(0.80, 0.042), fontsize=8,
                ha="center", va="top", arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6, shrinkB=4))
    ax.text(0.75, gw[3] + 0.003, f"+{gw[3]:.4f}", fontsize=8, ha="center", va="bottom", color=INK2)
    ax.text(0.0, -0.004, "base", fontsize=8, ha="left", va="top", color=INK2)
    ax.set_xlim(-0.06, 1.08)
    ax.set_ylim(-0.010, 0.095)
    ax.set_xticks(ws)
    ax.set_xticklabels(["0", ".25", ".5", ".75", "1"])
    ax.grid(axis="y")
    ax.set_xlabel("H1GW mixing weight w")
    ax.set_ylabel("composite gain over base (H1, real scorer)")
    ax.set_title("same-line data (H1GW)", fontsize=8.5)
    tag(ax, "a", x=-0.12, y=1.10)

    # ---- panel b: gains on one axis ------------------------------------------------------------------------------
    bx = fig.add_axes([0.575, 0.20, 0.41, 0.61])
    rows = ["one more source context", "M3 (control-PC deflation)", "M2 (agreement amplitude)",
            "best learned architecture ‡", "same-line H1GW (w = 1)"]
    y = {r: float(i) for i, r in enumerate(rows)}
    yA = y["best learned architecture ‡"]
    bx.axhspan(yA - 0.45, yA + 0.45, color=GRAY_TINT, lw=0, zorder=0)
    for r, v in [("M3 (control-PC deflation)", gM3), ("M2 (agreement amplitude)", gM2),
                 ("same-line H1GW (w = 1)", gw[-1])]:   # single H1 runs: emission-seed band
        bx.fill_between([v - spread / 2, v + spread / 2], y[r] - 0.22, y[r] + 0.22,
                        color=ORANGE_TINT if r.startswith("same") else BLUE_TINT, lw=0, zorder=1)
    ci_h(bx, y["one more source context"], pcv, pcv_ci[0], pcv_ci[1], S2, marker="^", ms=MS)
    bx.plot([gM3], [y["M3 (control-PC deflation)"]], marker="o", ms=MS, color=S1, mec="white", mew=0.8, ls="none",
            zorder=4)
    bx.plot([gM2], [y["M2 (agreement amplitude)"]], marker="o", ms=MS, color=S1, mec="white", mew=0.8, ls="none",
            zorder=4)
    ci_h(bx, yA, tl["mean"], tl["ci_lo"], tl["ci_hi"], S1, marker="D", ms=MS - 0.5, lw=1.4)
    bx.plot([sens], [yA], marker="D", ms=MS - 0.5, mfc="white", mec=S1, mew=1.1, ls="none", zorder=4)
    bx.plot([gw[-1]], [y["same-line H1GW (w = 1)"]], marker="o", ms=MS, color=S2, mec="white", mew=0.8, ls="none",
            zorder=4)
    lab = {"one more source context": (f"+{pcv:.5f} [{pcv_ci[0]:+.5f}, {pcv_ci[1]:+.5f}]", pcv_ci[1]),
           "M3 (control-PC deflation)": (f"+{gM3:.5f}", gM3 + spread / 2),
           "M2 (agreement amplitude)": (f"+{gM2:.4f}", gM2 + spread / 2)}
    for r, (t, x) in lab.items():
        bx.text(x + 0.0015, y[r], t, fontsize=8, va="center", ha="left", color=INK2)
    bx.text(tl["ci_lo"], yA + 0.27, f"seed mean +{tl['mean']:.4f} [{tl['ci_lo']:+.4f}, {tl['ci_hi']:+.4f}]",
            fontsize=8, va="center", ha="left", color=INK2)
    bx.text(sens + 0.0015, yA - 0.22, f"open: single line ({sens_line}) +{sens:.4f}", fontsize=8, va="center", ha="left",
            color=INK2)
    bx.text(gw[-1] + spread / 2 + 0.0015, y["same-line H1GW (w = 1)"], f"+{gw[-1]:.4f}", fontsize=8, va="center",
            ha="left", color=INK2)
    # ratio bracket
    yb = y["same-line H1GW (w = 1)"] - 0.36
    bx.annotate("", xy=(gw[-1], yb), xytext=(tl["mean"], yb),
                arrowprops=dict(arrowstyle="<|-|>", color=INK, lw=0.7, mutation_scale=6, shrinkA=0, shrinkB=0))
    bx.text(gw[-1] - 0.003, y["same-line H1GW (w = 1)"] + 0.02,
            f"≈ {ratio:.1f}× the largest modelling gain\n→ H4 confirmed", fontsize=8, ha="right", va="center",
            color=INK, fontweight="bold", linespacing=1.1)
    bx.axvline(0, color=INK2, lw=0.7, zorder=1)
    bx.set_yticks([y[r] for r in rows])
    bx.set_yticklabels(rows)
    bx.tick_params(axis="y", length=0)
    bx.spines["left"].set_visible(False)
    bx.grid(axis="x")
    bx.set_xlim(-0.004, 0.083)
    bx.set_xticks([0, 0.02, 0.04, 0.06, 0.08])
    bx.set_ylim(-0.6, len(rows) - 0.45)
    bx.set_xlabel("composite gain (linear axis)")
    bx.set_title("data value vs modelling gains", fontsize=8.5)
    tag(bx, "b", x=-0.60, y=1.10)
    h = [Line2D([], [], color=S2, marker="o", ms=MS, ls="none", label="data"),
         Line2D([], [], color=S1, marker="o", ms=MS, ls="none", label="modelling"),
         FancyBboxPatch((0, 0), 1, 1, color=BLUE_TINT, lw=0, label=f"band: single H1 run, emission-seed range ≈ {spread:.4f}")]
    fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.30, 0.995), ncol=3, fontsize=8, handlelength=1.3,
               columnspacing=1.2)
    fig.text(0.985, 0.015, "‡ K562/HepG2 architecture screen (own base, emission and seeds 0–2), not the H1 panel  ·  "
             "▲ offline LOFO, sum-side composite", fontsize=8, ha="right", va="bottom", color=INK2)
    save(fig, "fig6")


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, fig6):
        f()
