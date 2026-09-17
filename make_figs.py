import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
import os
os.makedirs("figures", exist_ok=True)
J = json.load(open("results/numbers.json"))
plt.rcParams.update({"font.family": "serif", "font.serif": ["Liberation Serif"], "mathtext.fontset": "stix",
                     "font.size": 8, "axes.linewidth": 0.6, "pdf.fonttype": 42, "savefig.bbox": "tight",
                     "savefig.pad_inches": 0.02})
COL = {"URAP": "#4E79A7", "EMAP": "#F28E2B", "SASI": "#59A14F", "M2AP": "#B07AA1", "LMAP": "#499894", "Ascon-Hash256": "#E15759"}
ORDER = sorted(J["qes"], key=lambda n: J["qes_rank"][n])
lab = lambda n: n.replace("-", "-\n")
def save(fig, s):
    fig.savefig("figures/" + s + ".pdf"); fig.savefig("figures/" + s + ".png", dpi=600); plt.close(fig)
W = 3.5; xs = np.arange(6)

# Fig 3 QES
fig, ax = plt.subplots(figsize=(W, 2.3)); q = [J["qes"][n] for n in ORDER]
ax.bar(xs, q, 0.62, color=[COL[n] for n in ORDER], edgecolor="black", lw=0.5, zorder=3)
for x, v in zip(xs, q): ax.text(x, v + 0.015, f"{v:.3f}", ha="center", fontsize=7)
for y in (0.2, 0.4, 0.6, 0.8): ax.axhline(y, color="0.6", lw=0.5, ls=(0, (4, 3)), zorder=1)
ax.set_ylim(0, 1.08); ax.set_xticks(xs); ax.set_xticklabels([lab(n) for n in ORDER], fontsize=7)
ax.set_ylabel("Quantum Exposure Score (QES)"); ax.set_xlabel("Authentication protocol"); ax.spines[["top", "right"]].set_visible(False)
save(fig, "Fig3_qes")

# Fig 4 QC
fig, ax = plt.subplots(figsize=(W, 2.3)); qc = [J["circuits"][n]["QC"] for n in ORDER]
ax.bar(xs, qc, 0.62, color=[COL[n] for n in ORDER], edgecolor="black", lw=0.5, zorder=3); ax.set_yscale("log"); ax.set_ylim(5e2, 3e7)
for x, v in zip(xs, qc): ax.text(x, v * 1.3, f"{int(v):,}", ha="center", fontsize=6.5)
ax.grid(axis="y", color="0.88", lw=0.5, zorder=0); ax.set_xticks(xs); ax.set_xticklabels([lab(n) for n in ORDER], fontsize=7)
ax.set_ylabel(r"Quantum cost, $QC = D \times W$"); ax.set_xlabel("Authentication protocol"); ax.spines[["top", "right"]].set_visible(False)
save(fig, "Fig4_qc")

# Fig 5 dual ranking
fig, ax = plt.subplots(figsize=(W, 3.1))
ax.plot([0.5, 6.5], [0.5, 6.5], color="0.45", lw=0.8, ls=(0, (5, 3)))
off = {"M2AP": (8, -2), "EMAP": (8, 3), "LMAP": (8, -2), "URAP": (8, 4), "SASI": (7, -16), "Ascon-Hash256": (-8, 4)}
for n in ORDER:
    x, y = J["qes_rank"][n], J["cost_rank"][n]
    if x != y: ax.plot([x, x], [x, y], color=COL[n], lw=0.8, ls=":")
    ax.scatter(x, y, s=46, color=COL[n], edgecolor="black", lw=0.5, zorder=4)
    dx, dy = off[n]
    ax.annotate(f"{n}\nQES {J['qes'][n]:.3f}, $2^{{{J['attack'][n]['log2_Ttotal']:.1f}}}$", (x, y), xytext=(dx, dy),
                textcoords="offset points", ha="left" if dx > 0 else "right", va="center", fontsize=6.4, linespacing=1.15)
ax.text(0.7, 6.25, rf"Kendall $\tau = {J['tau']:.3f}$", fontsize=7.5)
ax.set_xlim(0.5, 6.5); ax.set_ylim(0.5, 6.5); ax.set_xticks(range(1, 7)); ax.set_yticks(range(1, 7)); ax.set_aspect("equal")
ax.set_xlabel("Oracle-construction rank (QES; 1 = most exposed)"); ax.set_ylabel(r"Attack-cost rank ($T_{\mathrm{total}}$; 1 = cheapest)")
ax.grid(color="0.9", lw=0.5); save(fig, "Fig5_dual_ranking")

# Fig 6 sensitivity heatmap
dist = np.array([J["sens"]["dist"][n] for n in ORDER])
fig, ax = plt.subplots(figsize=(W, 2.45)); im = ax.imshow(dist, cmap="Blues", vmin=0, vmax=100, aspect="auto")
for i in range(6):
    for j in range(6):
        v = dist[i, j]; t = "\u2013" if v == 0 else ("<0.1" if v < 0.1 else f"{v:.1f}")
        ax.text(j, i, t, ha="center", va="center", fontsize=6.5, color="white" if v > 55 else "black")
    ax.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False, ec="black", lw=1.2, clip_on=False))
ax.set_xticks(range(6)); ax.set_xticklabels(range(1, 7)); ax.set_yticks(range(6))
ax.set_yticklabels([f"{n} ({i+1})" for i, n in enumerate(ORDER)]); ax.tick_params(length=0)
ax.set_xlabel("Rank position under random weights"); ax.set_ylabel("Protocol (equal-weight rank)")
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.set_label("Share of 50,000 weightings (%)", fontsize=7); cb.ax.tick_params(labelsize=6.5)
save(fig, "Fig8_rank_stability")

# Fig 7 oracle budget frontier (both depth assumptions)
fig, ax = plt.subplots(figsize=(W, 2.7))
off = {"URAP": (6, 0), "EMAP": (6, -4), "SASI": (6, 2), "M2AP": (6, -2), "LMAP": (-6, 5), "Ascon-Hash256": (-6, 5)}
for n in ORDER:
    w = J["attack"][n]["W_oracle"]; dp, ds = J["attack"][n]["depth_iter"], J["attack_serial"][n]["depth_iter"]
    if ds != dp: ax.plot([w, w], [dp, ds], color=COL[n], lw=0.9, zorder=2)
    ax.scatter(w, dp, s=38, facecolor="white", edgecolor=COL[n], lw=1.3, zorder=3)
    ax.scatter(w, ds, s=38, color=COL[n], edgecolor="black", lw=0.5, zorder=4)
    dx, dy = off[n]
    ax.annotate(n, (w, np.sqrt(dp * ds)), xytext=(dx, dy), textcoords="offset points", fontsize=6.6, ha="left" if dx > 0 else "right", va="center")
ax.scatter([], [], s=30, color="0.4", edgecolor="black", lw=0.5, label="All operations in series")
ax.scatter([], [], s=30, facecolor="white", edgecolor="0.4", lw=1.2, label="Messages in parallel")
ax.legend(frameon=False, fontsize=6.5, loc="upper right")
ax.set_xscale("log"); ax.set_yscale("log")
ax.xaxis.set_major_locator(FixedLocator([300, 500, 1000, 1500])); ax.xaxis.set_major_formatter(ScalarFormatter()); ax.xaxis.set_minor_formatter(NullFormatter())
ax.yaxis.set_major_locator(FixedLocator([500, 1000, 2000, 5000, 10000])); ax.yaxis.set_major_formatter(ScalarFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_xlim(280, 1700); ax.set_ylim(500, 20000); ax.grid(which="both", color="0.92", lw=0.5)
ax.set_xlabel(r"Oracle width $W_{\mathrm{oracle}}$ (logical qubits)"); ax.set_ylabel(r"$T$-depth per Grover iteration")
save(fig, "Fig6_budget_frontier")

# Fig 8 NIST margin vs word size
fig, axes = plt.subplots(1, 3, figsize=(7.16, 2.3), sharey=True)
ulw = [n for n in ORDER if n != "Ascon-Hash256"]
for ax, md in zip(axes, ("40", "64", "96")):
    for n in ulw:
        bs = sorted(map(int, J["sweep"][n])); ys = [J["sweep"][n][str(b)][md] for b in bs]
        ax.plot(bs, ys, marker="o", ms=2.5, lw=1.1, color=COL[n], label=n)
    ax.axhline(0, color="black", lw=0.8); ax.axvline(32, color="0.6", lw=0.6, ls=":")
    ax.set_ylim(-110, 160); ax.set_xticks([16, 32, 64, 96, 128]); ax.grid(color="0.9", lw=0.5)
    ax.set_title(f"MAXDEPTH $2^{{{md}}}$", fontsize=8); ax.set_xlabel("Word size $b$ (bits)")
axes[0].set_ylabel("Margin to NIST category 1 (bits)")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=5, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.2))
save(fig, "Fig7_word_size")
print("ok", ORDER)
