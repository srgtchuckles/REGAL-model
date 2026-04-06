"""
GPS (Galinpepimut-S) · REGAL · Mixture Cure Model
===================================================
Parametric survival analysis with HLA-stratified cure fraction estimation.
Grounded in Phase 2 immunologic correlate data (Blood 2019) and historical
AML CR2 BAT benchmarks.

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.special import erf
from dataclasses import dataclass


# ─────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────

@dataclass
class ModelParams:
    # HLA / Immunologic parameters
    hla_positive_rate: float = 0.45        # ~45% HLA-A*02:01+ in European AML population
    cd8_response_rate: float = 0.86        # 6/7 HLA-A*02+ tested → CD8+ response (Phase 2)
    cd4_response_rate: float = 0.44        # 4/9 patients → CD4+ response (Phase 2)
    hla_neg_pop: float = 0.55              # 
    cure_fraction_responders: float = 0.575 # ~57.5% DFS plateau in immunologic responders (Fig 5)
    cure_fraction_hla_neg: float = 0.3
    hla_neg_cure_fraction: float = 0.025   # Small residual cure via non-A*02 MHC I + CD4 epitopes

    # BAT arm — historical AML CR2 (not eligible for transplant)
    bat_median_os: float = 11.0            # months; literature range 8–16mo
    bat_weibull_shape: float = 1.2         # k > 1 = increasing hazard (appropriate for AML)

    # GPS arm — non-cured susceptible subpopulation
    gps_median_os: float = 20.0            # months for non-cured GPS patients
    gps_weibull_shape: float = 1.3

    # Trial
    n_patients: int = 116
    follow_up_months: int = 72


# ─────────────────────────────────────────────
# Survival Functions
# ─────────────────────────────────────────────

def weibull_lambda(median_os: float, shape: float) -> float:
    """Scale parameter from median and shape."""
    return median_os / (np.log(2) ** (1 / shape))


def weibull_survival(t: np.ndarray, median_os: float, shape: float) -> np.ndarray:
    """S(t) = exp(-(t/lambda)^k)"""
    lam = weibull_lambda(median_os, shape)
    return np.exp(-((t / lam) ** shape))


def mixture_cure_survival(
    t: np.ndarray,
    cure_fraction: float,
    median_os: float,
    shape: float
) -> np.ndarray:
    """
    Mixture cure model:
        S(t) = pi + (1 - pi) * S_susceptible(t)
    where pi = cure fraction, S_susceptible = Weibull for non-cured patients.
    As t → inf, S(t) → pi (the plateau).
    """
    s_susceptible = weibull_survival(t, median_os, shape)
    return cure_fraction + (1 - cure_fraction) * s_susceptible


def compute_aggregate_cure_fraction(p: ModelParams) -> float:
    """
    Derive aggregate cure fraction from HLA-stratified components.

    Chain:
      HLA-A*02:01+ × CD8 response rate × DFS plateau (cure|responder)
      + HLA-A*02:01+ × CD8 resp × CD4 bonus (durability amplifier)
      + HLA-A*02:01- × residual cure (non-A*02 epitopes)
    """
    cd8_cure = p.hla_positive_rate * p.cd8_response_rate * p.cure_fraction_responders
    hla_negative_cure = p.cd4_response_rate * p.hla_neg_pop * p.cure_fraction_hla_neg
    return min(cd8_cure + hla_negative_cure, 0.45)


# ─────────────────────────────────────────────
# Hazard Rate (instantaneous)
# ─────────────────────────────────────────────

def weibull_hazard(t: np.ndarray, median_os: float, shape: float) -> np.ndarray:
    """h(t) = (k/lambda) * (t/lambda)^(k-1)"""
    lam = weibull_lambda(median_os, shape)
    t_safe = np.maximum(t, 1e-9)
    return (shape / lam) * ((t_safe / lam) ** (shape - 1))


def mixture_cure_hazard(
    t: np.ndarray,
    cure_fraction: float,
    median_os: float,
    shape: float
) -> np.ndarray:
    """
    h_mix(t) = [(1 - pi) * f_susceptible(t)] / S_mix(t)
    where f = h * S for the susceptible subpopulation.
    """
    s_sus = weibull_survival(t, median_os, shape)
    h_sus = weibull_hazard(t, median_os, shape)
    f_sus = h_sus * s_sus  # density for susceptible
    s_mix = mixture_cure_survival(t, cure_fraction, median_os, shape)
    return np.where(s_mix > 1e-9, (1 - cure_fraction) * f_sus / s_mix, 0.0)


# ─────────────────────────────────────────────
# Event Velocity (6-month windows)
# ─────────────────────────────────────────────

def event_velocity(
    t: np.ndarray,
    survival_fn,
    n_patients: int,
    window: int = 6
) -> tuple:
    """
    Expected events in each [t, t+window] interval.
    Events = N * [S(t) - S(t+window)]
    Returns (window_starts, event_counts).
    """
    starts = np.arange(0, t[-1] - window + 1, window)
    counts = []
    for s in starts:
        s1 = survival_fn(np.array([s]))[0]
        s2 = survival_fn(np.array([s + window]))[0]
        counts.append(n_patients * (s1 - s2))
    return starts, np.array(counts)


# ─────────────────────────────────────────────
# Approximate Log-Rank Statistics
# ─────────────────────────────────────────────

def approximate_hazard_ratio(
    p: ModelParams,
    cure_fraction: float,
    eval_time: float = None
) -> float:
    """
    Approximate HR at the BAT median OS timepoint.
    HR = log(S_GPS) / log(S_BAT) — ratio of cumulative hazards.
    """
    t = eval_time or p.bat_median_os
    s_bat = weibull_survival(np.array([t]), p.bat_median_os, p.bat_weibull_shape)[0]
    s_gps = mixture_cure_survival(
        np.array([t]), cure_fraction, p.gps_median_os, p.gps_weibull_shape
    )[0]
    if s_bat <= 0 or s_gps <= 0:
        return np.nan
    return np.log(s_gps) / np.log(s_bat)


def normal_cdf(x: float) -> float:
    return 0.5 * (1 + erf(x / np.sqrt(2)))


def approximate_pvalue(hr: float, n_events: int) -> float:
    """Log-rank z-score approximation."""
    z = -np.log(hr) * np.sqrt(n_events / 4)
    return min(1 - normal_cdf(z), 0.5)


# ─────────────────────────────────────────────
# Sensitivity Analysis
# ─────────────────────────────────────────────

def sensitivity_analysis(base_params: ModelParams):
    """
    Sweep BAT median OS (8–16mo) and cure fraction (10–25%).
    Returns grid of HRs and p-values.
    """
    bat_range = np.linspace(8, 16, 9)
    cure_range = np.linspace(0.10, 0.25, 9)

    hr_grid = np.zeros((len(cure_range), len(bat_range)))
    pval_grid = np.zeros_like(hr_grid)
    n_events = int(base_params.n_patients * 0.65)

    for i, cf in enumerate(cure_range):
        for j, bat_mos in enumerate(bat_range):
            p = ModelParams(**{**base_params.__dict__,
                               'bat_median_os': bat_mos})
            t_eval = bat_mos
            s_bat = weibull_survival(
                np.array([t_eval]), bat_mos, p.bat_weibull_shape)[0]
            s_gps = mixture_cure_survival(
                np.array([t_eval]), cf, p.gps_median_os, p.gps_weibull_shape)[0]
            hr = np.log(s_gps) / np.log(s_bat) if s_bat > 0 and s_gps > 0 else np.nan
            hr_grid[i, j] = hr
            pval_grid[i, j] = approximate_pvalue(hr, n_events) if not np.isnan(hr) else np.nan

    return bat_range, cure_range, hr_grid, pval_grid


# ─────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────

DARK_BG   = "#060f1e"
PANEL_BG  = "#0a1628"
GRID_COL  = "#0f2035"
BLUE      = "#3b82f6"
ORANGE    = "#f97316"
HIST_LOW  = "#1a3050"
HIST_HIGH = "#2a4a6a"
TEXT_DIM  = "#5a7a9a"
TEXT_MID  = "#8aaac8"
TEXT_HI   = "#e2eaf5"
GREEN     = "#22c55e"
AMBER     = "#f59e0b"
RED       = "#ef4444"


def style_ax(ax, title=None):
    ax.set_facecolor(PANEL_BG)
    ax.tick_params(colors=TEXT_DIM, labelsize=8)
    ax.xaxis.label.set_color(TEXT_DIM)
    ax.yaxis.label.set_color(TEXT_DIM)
    for spine in ax.spines.values():
        spine.set_edgecolor(GRID_COL)
    ax.grid(True, color=GRID_COL, linewidth=0.5, linestyle="--", alpha=0.7)
    if title:
        ax.set_title(title, color=TEXT_MID, fontsize=9,
                     fontfamily="monospace", pad=10, loc="left")


def plot_model(p: ModelParams = None):
    if p is None:
        p = ModelParams()

    t = np.linspace(0, p.follow_up_months, 500)
    cure_fraction = compute_aggregate_cure_fraction(p)
    hr = approximate_hazard_ratio(p, cure_fraction)
    n_events = int(p.n_patients * 0.65)
    pval = approximate_pvalue(hr, n_events)

    # Survival curves
    s_bat = weibull_survival(t, p.bat_median_os, p.bat_weibull_shape)
    s_gps = mixture_cure_survival(t, cure_fraction, p.gps_median_os, p.gps_weibull_shape)
    s_hist_low = weibull_survival(t, 8, 1.1)
    s_hist_high = weibull_survival(t, 16, 1.3)

    # Hazard curves
    h_bat = weibull_hazard(t[1:], p.bat_median_os, p.bat_weibull_shape)
    h_gps = mixture_cure_hazard(t[1:], cure_fraction, p.gps_median_os, p.gps_weibull_shape)

    # Event velocity
    vel_bat_starts, vel_bat = event_velocity(
        t, lambda x: weibull_survival(x, p.bat_median_os, p.bat_weibull_shape), p.n_patients)
    vel_gps_starts, vel_gps = event_velocity(
        t, lambda x: mixture_cure_survival(x, cure_fraction, p.gps_median_os, p.gps_weibull_shape), p.n_patients)

    # Sensitivity
    bat_range, cure_range, hr_grid, pval_grid = sensitivity_analysis(p)

    # ── Figure layout ──
    fig = plt.figure(figsize=(18, 13), facecolor=DARK_BG)
    fig.suptitle(
        "GPS · REGAL  |  Mixture Cure Model",
        color=TEXT_HI, fontsize=16, fontfamily="monospace", y=0.98, x=0.02, ha="left"
    )
    fig.text(0.02, 0.955, "Galinpepimut-S · AML CR2 · Phase 3",
             color=TEXT_DIM, fontsize=9, fontfamily="monospace")

    gs = gridspec.GridSpec(3, 3, figure=fig,
                           hspace=0.45, wspace=0.35,
                           left=0.06, right=0.97, top=0.93, bottom=0.06)

    # ── Stat banner ──
    # fix spacing between the variables here
    stats = [
        ("Cure Fraction π", f"{cure_fraction*100:.1f}%", BLUE),
        ("Hazard Ratio",    f"{hr:.3f}",                  GREEN if hr < 0.7 else AMBER),
        ("Est. p-value",    f"{'<0.001' if pval < 0.001 else f'{pval:.3f}'}",
                            GREEN if pval < 0.05 else RED),
        ("GPS 36mo OS",     f"{mixture_cure_survival(np.array([36]), cure_fraction, p.gps_median_os, p.gps_weibull_shape)[0]*100:.1f}%", TEXT_MID),
        ("BAT 36mo OS",     f"{weibull_survival(np.array([36]), p.bat_median_os, p.bat_weibull_shape)[0]*100:.1f}%", ORANGE),
        ("GPS 60mo OS",     f"{mixture_cure_survival(np.array([60]), cure_fraction, p.gps_median_os, p.gps_weibull_shape)[0]*100:.1f}%", TEXT_MID),
    ]
    for idx, (lbl, val, col) in enumerate(stats):
        x = 0.1 + idx * 0.07
        fig.text(x, 0.925, lbl, color=TEXT_DIM, fontsize=7, fontfamily="monospace", transform=fig.transFigure)
        fig.text(x, 0.908, val, color=col,      fontsize=13, fontfamily="monospace",
                 fontweight="bold", transform=fig.transFigure)

    # ── Panel 1: Survival curves ──
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.fill_between(t, s_hist_low * 100, s_hist_high * 100,
                     color=HIST_LOW, alpha=0.4, label="BAT historical range (8–16mo)")
    ax1.plot(t, s_bat * 100,  color=ORANGE, lw=2.5, label=f"BAT modeled (mOS={p.bat_median_os}mo)")
    ax1.plot(t, s_gps * 100,  color=BLUE,   lw=3.0, label=f"GPS mixture cure (π={cure_fraction*100:.1f}%)")
    ax1.axhline(cure_fraction * 100, color=BLUE, lw=1, ls=":", alpha=0.6,
                label=f"Cure plateau π={cure_fraction*100:.1f}%")
    ax1.set_xlabel("Months"); ax1.set_ylabel("Overall Survival (%)")
    ax1.set_ylim(0, 105); ax1.set_xlim(0, p.follow_up_months)
    ax1.legend(fontsize=7, facecolor=DARK_BG, edgecolor=GRID_COL,
               labelcolor=TEXT_MID, loc="upper right")
    style_ax(ax1, "OVERALL SURVIVAL  ·  Mixture Cure vs BAT")

    # ── Panel 2: Hazard rates ──
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.plot(t[1:], h_bat, color=ORANGE, lw=2, label="BAT")
    ax2.plot(t[1:], h_gps, color=BLUE,   lw=2, label="GPS")
    ax2.set_xlabel("Months"); ax2.set_ylabel("Hazard Rate h(t)")
    ax2.set_xlim(0, p.follow_up_months)
    ax2.legend(fontsize=7, facecolor=DARK_BG, edgecolor=GRID_COL, labelcolor=TEXT_MID)
    style_ax(ax2, "INSTANTANEOUS HAZARD")

    # ── Panel 3: Event velocity ──
    ax3 = fig.add_subplot(gs[1, :2])
    width = 2.2
    labels = [f"Mo {int(s)+1}–{int(s)+6}" for s in vel_bat_starts]
    x_pos = np.arange(len(labels))
    ax3.bar(x_pos - width/2, vel_bat, width, color=ORANGE, alpha=0.85, label="BAT")
    ax3.bar(x_pos + width/2, vel_gps, width, color=BLUE,   alpha=0.85, label="GPS")
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
    ax3.set_ylabel("Expected Events")
    ax3.legend(fontsize=7, facecolor=DARK_BG, edgecolor=GRID_COL, labelcolor=TEXT_MID)
    style_ax(ax3, "EVENT VELOCITY  ·  6-Month Windows")

    # ── Panel 4: Cure fraction derivation ──
    ax4 = fig.add_subplot(gs[1, 2])
    ax4.set_facecolor(PANEL_BG)
    ax4.axis("off")
    cd8_comp = p.hla_positive_rate * p.cd8_response_rate * p.cure_fraction_responders
    hla_neg_comp = p.cd4_response_rate * p.hla_neg_pop * p.cure_fraction_hla_neg
    lines = [
        ("CURE FRACTION DERIVATION", TEXT_MID, 10, True),
        ("", TEXT_DIM, 8, False),
        (f"HLA-A*02:01+ rate:  {p.hla_positive_rate*100:.0f}%", TEXT_DIM, 8, False),
        (f"CD8+ resp rate:     {p.cd8_response_rate*100:.0f}%", TEXT_DIM, 8, False),
        (f"Cure | responder:   {p.cure_fraction_responders*100:.0f}%", TEXT_DIM, 8, False),
        ("", TEXT_DIM, 8, False),
        (f"HLA-A*02:01- rate: {p.hla_neg_pop*100:.0f}%", TEXT_DIM, 8, False),
        (f"CD4+ resp rate: {p.cd4_response_rate*100:.0f}%", TEXT_DIM, 8, False),
        (f"Cure | long pep resp: {p.cure_fraction_hla_neg*100:.0f}%", TEXT_DIM, 8, False),
        ("", TEXT_DIM, 8, False),
        (f"CD8 component:  {cd8_comp*100:.1f}%", BLUE, 9, False),
        (f"HLA neg component: {hla_neg_comp*100:.1f}%", BLUE, 9, False),
        ("─" * 22, TEXT_DIM, 8, False),
        (f"Aggregate π:    {cure_fraction*100:.1f}%", BLUE, 12, True),
        ("", TEXT_DIM, 8, False),
        (f"HR (approx):    {hr:.3f}", GREEN if hr < 0.7 else AMBER, 9, False),
        (f"p-value:        {'<0.001' if pval < 0.001 else f'{pval:.3f}'}", GREEN if pval < 0.05 else RED, 9, False),
    ]
    y = 0.95
    for text, color, size, bold in lines:
        ax4.text(0.05, y, text, color=color, fontsize=size,
                 fontfamily="monospace", fontweight="bold" if bold else "normal",
                 transform=ax4.transAxes)
        y -= 0.072
    for spine in ax4.spines.values():
        spine.set_edgecolor(GRID_COL)
    style_ax(ax4)

    # ── Panel 5: HR sensitivity heatmap (BAT mOS vs cure fraction) ──
    ax5 = fig.add_subplot(gs[2, :2])
    im = ax5.imshow(hr_grid, aspect="auto", origin="lower",
                    cmap="RdYlGn_r", vmin=0.4, vmax=1.1,
                    extent=[bat_range[0], bat_range[-1],
                            cure_range[0] * 100, cure_range[-1] * 100])
    ax5.set_xlabel("BAT Median OS (months)")
    ax5.set_ylabel("Cure Fraction π (%)")

    # Significance contour at HR=0.75
    cs = ax5.contour(bat_range, cure_range * 100, hr_grid,
                     levels=[0.65, 0.75, 0.85],
                     colors=["white", "lightyellow", "orange"],
                     linewidths=[1.5, 1.0, 0.8], linestyles=["--", ":", ":"])
    ax5.clabel(cs, fmt="HR=%.2f", fontsize=7, colors="white")

    # Mark current params
    ax5.scatter([p.bat_median_os], [cure_fraction * 100],
                color="white", s=80, zorder=5, marker="x", linewidths=2)
    ax5.text(p.bat_median_os + 0.2, cure_fraction * 100 + 0.3,
             "current", color="white", fontsize=7, fontfamily="monospace")

    plt.colorbar(im, ax=ax5, label="Hazard Ratio", shrink=0.8).ax.yaxis.label.set_color(TEXT_DIM)
    style_ax(ax5, "SENSITIVITY  ·  HR across BAT mOS × Cure Fraction  (✕ = current params)")

    # ── Panel 6: p-value sensitivity ──
    ax6 = fig.add_subplot(gs[2, 2])
    im2 = ax6.imshow(pval_grid, aspect="auto", origin="lower",
                     cmap="RdYlGn_r", vmin=0, vmax=0.2,
                     extent=[bat_range[0], bat_range[-1],
                             cure_range[0] * 100, cure_range[-1] * 100])
    ax6.set_xlabel("BAT Median OS (months)")
    ax6.set_ylabel("Cure Fraction π (%)")
    cs2 = ax6.contour(bat_range, cure_range * 100, pval_grid,
                      levels=[0.05], colors=["white"], linewidths=[1.5], linestyles=["--"])
    ax6.clabel(cs2, fmt="p=%.2f", fontsize=7, colors="white")
    ax6.scatter([p.bat_median_os], [cure_fraction * 100],
                color="white", s=80, zorder=5, marker="x", linewidths=2)
    plt.colorbar(im2, ax=ax6, label="p-value", shrink=0.8).ax.yaxis.label.set_color(TEXT_DIM)
    style_ax(ax6, "SENSITIVITY  ·  p-value")

    import os
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gps_model_output.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=DARK_BG)
    print(f"Saved: {out_path}")
    plt.show()
    return fig


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

if __name__ == "__main__":
    p = ModelParams()
    cure_fraction = compute_aggregate_cure_fraction(p)

    print("=" * 52)
    print("  GPS · REGAL · Mixture Cure Model")
    print("=" * 52)
    print(f"  HLA-A*02:01+ rate:         {p.hla_positive_rate*100:.0f}%")
    print(f"  CD8+ response rate:         {p.cd8_response_rate*100:.0f}%")
    print(f"  CD4+ response rate:         {p.cd4_response_rate*100:.0f}%")
    print(f"  Cure | responder (Ph2 KM):  {p.cure_fraction_responders*100:.0f}%")
    print(f"  HLA- residual cure:         {p.hla_neg_cure_fraction*100:.0f}%")
    print(f"  BAT median OS:              {p.bat_median_os}mo")
    print(f"  GPS median OS (non-cured):  {p.gps_median_os}mo")
    print("-" * 52)
    print(f"  Aggregate cure fraction π:  {cure_fraction*100:.1f}%")

    hr = approximate_hazard_ratio(p, cure_fraction)
    n_events = int(p.n_patients * 0.65)
    pval = approximate_pvalue(hr, n_events)

    print(f"  Approx Hazard Ratio:        {hr:.3f}")
    print(f"  Approx p-value:             {'<0.001' if pval < 0.001 else f'{pval:.4f}'}")

    t_pts = [12, 24, 36, 48, 60]
    print("\n  Survival estimates:")
    print(f"  {'Month':<8} {'GPS':>8} {'BAT':>8} {'Δ':>8}")
    print(f"  {'─'*36}")
    for t in t_pts:
        s_g = mixture_cure_survival(np.array([t]), cure_fraction, p.gps_median_os, p.gps_weibull_shape)[0]
        s_b = weibull_survival(np.array([t]), p.bat_median_os, p.bat_weibull_shape)[0]
        print(f"  {t:<8} {s_g*100:>7.1f}% {s_b*100:>7.1f}% {(s_g-s_b)*100:>+7.1f}%")

    print("=" * 52)
    plot_model(p)
