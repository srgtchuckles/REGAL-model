"""
GPS · REGAL · Monte Carlo Enrollment Simulation
================================================
Simulates realistic trial enrollment over calendar time, draws
patient-level event times from the mixture cure model, and validates
the model's cure fraction estimate against the two SEC-filed event
count anchors:
    - 60 events at calendar month 46
    - 72 events at calendar month 58

The key distinction from the analytical model:
    - Analytical model: all patients start at t=0 simultaneously
    - This model: patients enroll over ~30 months (Jan 2020 – Jun 2022)
      and are followed until their event or the calendar cutoff

Each simulation run generates a synthetic trial. Running 10,000 runs
gives a distribution of expected event counts at each calendar anchor,
which we compare to the observed SEC-filed counts to assess whether
our cure fraction estimate is consistent with reality.

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from dataclasses import dataclass
from typing import Tuple
import os
from scipy.interpolate import interp1d


# ── import from main model ──────────────────────────────────────────
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from REGAL_model import (
    ModelParams,
    weibull_lambda,
    weibull_survival,
    mixture_cure_survival,
    compute_aggregate_cure_fraction,
)

# ── cosmetics (match main model) ────────────────────────────────────
DARK_BG   = "#060f1e"
PANEL_BG  = "#0a1628"
GRID_COL  = "#0f2035"
BLUE      = "#3b82f6"
ORANGE    = "#f97316"
TEXT_DIM  = "#5a7a9a"
TEXT_MID  = "#8aaac8"
TEXT_HI   = "#e2eaf5"
GREEN     = "#22c55e"
AMBER     = "#f59e0b"
RED       = "#ef4444"
PURPLE    = "#a855f7"


# ─────────────────────────────────────────────
# Trial Timeline Parameters
# ─────────────────────────────────────────────

@dataclass
class TrialTimeline:
    # Calendar months from trial start (Jan 2020 = month 0)
    enrollment_start: float = 0.0      # Jan 2020
    enrollment_end: float   = 30.0     # ~Jun 2022 (estimated ~2.5yr accrual)
    calendar_cutoff: float  = 72.0     # ~Jan 2026 (6yr from start)

    # SEC-filed event count anchors (calendar months from trial start)
    anchor_1_month:  float = 46.0      # Nov 2023
    anchor_1_events: int   = 60        # observed deaths
    anchor_2_month:  float = 58.0      # Nov 2024
    anchor_2_events: int   = 72        # observed deaths

    # Final analysis trigger
    final_analysis_events: int = 80

    # Randomization ratio
    gps_fraction: float = 0.50         # 1:1 randomization


# ─────────────────────────────────────────────
# Patient-Level Event Time Sampling
# ─────────────────────────────────────────────

def sample_weibull_time(median_os: float, shape: float, rng: np.random.Generator) -> float:
    """
    Sample event time from Weibull distribution via inverse CDF.
    If U ~ Uniform(0,1), then T = lambda * (-log(U))^(1/k)
    """
    lam = weibull_lambda(median_os, shape)
    u = rng.uniform()
    return lam * (-np.log(u + 1e-12)) ** (1.0 / shape)


def sample_patient_event_time(
    arm: str,
    cure_fraction: float,
    p: ModelParams,
    rng: np.random.Generator
) -> float:
    """
    Sample a single patient's time-to-event from enrollment.

    Returns:
        float: months from enrollment to event (inf = cured / no event)
    """
    if arm == "GPS":
        # First determine if patient is in the cured fraction
        if rng.uniform() < cure_fraction:
            return np.inf  # cured — will never have the event
        else:
            return sample_weibull_time(p.gps_median_os, p.gps_weibull_shape, rng)
    else:  # BAT — no cure fraction
        return sample_weibull_time(p.bat_median_os, p.bat_weibull_shape, rng)


# ─────────────────────────────────────────────
# Enrollment Distribution
# ─────────────────────────────────────────────

def sample_enrollment_times(
    n_patients: int,
    timeline: TrialTimeline,
    rng: np.random.Generator,
    pattern: str = "uniform"
) -> np.ndarray:
    """
    Sample enrollment times (months from trial start) for all patients.

    pattern options:
        "uniform"  — constant enrollment rate throughout accrual period
        "ramp"     — slow start, ramps up (common in real trials)
        "bell"     — peaks in middle of accrual
    """
    duration = timeline.enrollment_end - timeline.enrollment_start

    if pattern == "uniform":
        times = rng.uniform(0, duration, n_patients)

    elif pattern == "ramp":
        # Triangular distribution — slow start, faster later
        u = rng.uniform(0, 1, n_patients)
        times = duration * np.sqrt(u)  # inverse CDF of linear ramp

    elif pattern == "bell":
        # Normal centered at midpoint, truncated
        mid = duration / 2
        sd = duration / 4
        times = rng.normal(mid, sd, n_patients)
        times = np.clip(times, 0, duration)

    return np.sort(times + timeline.enrollment_start)


# ─────────────────────────────────────────────
# Single Trial Simulation
# ─────────────────────────────────────────────

def simulate_trial(
    p: ModelParams,
    timeline: TrialTimeline,
    cure_fraction: float,
    rng: np.random.Generator,
    enrollment_pattern: str = "uniform"
) -> dict:
    """
    Simulate a single REGAL trial run.

    Returns dict with:
        - enrollment_times: array of calendar enrollment months
        - arms: array of "GPS" or "BAT" per patient
        - event_times_from_enrollment: array of event times (inf = censored/cured)
        - calendar_event_times: enrollment + event_time
        - events_by_calendar_month: cumulative event count at each calendar month
        - n_events_at_anchor_1: events at month 46
        - n_events_at_anchor_2: events at month 58
        - month_of_80th_event: calendar month when 80th event occurs
    """
    n = p.n_patients
    n_gps = int(n * timeline.gps_fraction)
    n_bat = n - n_gps

    # Enrollment times
    enroll_times = sample_enrollment_times(n, timeline, rng, enrollment_pattern)

    # Arm assignment (first n_gps get GPS, rest get BAT — random order)
    arms = np.array(["GPS"] * n_gps + ["BAT"] * n_bat)
    rng.shuffle(arms)

    # Sample event times from enrollment for each patient
    event_times_from_enroll = np.array([
        sample_patient_event_time(arm, cure_fraction, p, rng)
        for arm in arms
    ])

    # Convert to calendar time
    calendar_event_times = enroll_times + event_times_from_enroll

    # Observed event times: min(event, cutoff)
    # Patient is censored if event > calendar cutoff OR event = inf
    calendar_cutoff = timeline.calendar_cutoff
    observed_calendar = np.minimum(calendar_event_times, calendar_cutoff)
    had_event = (calendar_event_times <= calendar_cutoff) & np.isfinite(calendar_event_times)

    # Cumulative events over calendar months
    event_calendar_times = np.sort(calendar_event_times[had_event & np.isfinite(calendar_event_times)])

    def events_by_month(cal_month):
        return int(np.sum(event_calendar_times <= cal_month))

    n_at_anchor_1 = events_by_month(timeline.anchor_1_month)
    n_at_anchor_2 = events_by_month(timeline.anchor_2_month)

    # Month of 80th event
    if len(event_calendar_times) >= timeline.final_analysis_events:
        month_80 = event_calendar_times[timeline.final_analysis_events - 1]
    else:
        month_80 = np.inf  # never reached in this sim

    # Follow-up time from enrollment (for KM reconstruction)
    follow_up = np.minimum(
        event_times_from_enroll,
        calendar_cutoff - enroll_times
    )
    follow_up = np.maximum(follow_up, 0)

    return {
        "enrollment_times":         enroll_times,
        "arms":                     arms,
        "event_times_from_enroll":  event_times_from_enroll,
        "calendar_event_times":     calendar_event_times,
        "had_event":                had_event,
        "follow_up":                follow_up,
        "n_at_anchor_1":            n_at_anchor_1,
        "n_at_anchor_2":            n_at_anchor_2,
        "month_of_80th_event":      month_80,
        "n_gps":                    n_gps,
        "n_bat":                    n_bat,
    }


# ─────────────────────────────────────────────
# Monte Carlo Runner
# ─────────────────────────────────────────────

def run_monte_carlo(
    p: ModelParams,
    timeline: TrialTimeline,
    cure_fraction: float,
    n_sims: int = 10_000,
    seed: int = 42,
    enrollment_pattern: str = "uniform"
) -> dict:
    """
    Run n_sims trial simulations and collect summary statistics.
    """
    rng = np.random.default_rng(seed)

    anchor_1_counts = np.zeros(n_sims, dtype=int)
    anchor_2_counts = np.zeros(n_sims, dtype=int)
    month_80_vals   = np.zeros(n_sims)

    print(f"  Running {n_sims:,} simulations...", end="", flush=True)

    for i in range(n_sims):
        result = simulate_trial(p, timeline, cure_fraction, rng, enrollment_pattern)
        anchor_1_counts[i] = result["n_at_anchor_1"]
        anchor_2_counts[i] = result["n_at_anchor_2"]
        month_80_vals[i]   = result["month_of_80th_event"]

        if (i + 1) % 2000 == 0:
            print(f" {i+1//1000}k", end="", flush=True)

    print(" done.")

    # p-value: fraction of sims with events >= observed (one-sided upper)
    # We want to know if observed counts are consistent with our model
    p_anchor_1 = np.mean(anchor_1_counts >= timeline.anchor_1_events)
    p_anchor_2 = np.mean(anchor_2_counts >= timeline.anchor_2_events)

    # Two-sided consistency: is observed within central 95% of simulated?
    ci_1 = np.percentile(anchor_1_counts, [2.5, 97.5])
    ci_2 = np.percentile(anchor_2_counts, [2.5, 97.5])

    # Filter inf for month_80 stats
    finite_80 = month_80_vals[np.isfinite(month_80_vals)]

    return {
        "anchor_1_counts":  anchor_1_counts,
        "anchor_2_counts":  anchor_2_counts,
        "month_80_vals":    month_80_vals,
        "finite_80":        finite_80,
        # Anchor 1
        "a1_mean":  anchor_1_counts.mean(),
        "a1_median":np.median(anchor_1_counts),
        "a1_ci":    ci_1,
        "a1_p":     p_anchor_1,
        # Anchor 2
        "a2_mean":  anchor_2_counts.mean(),
        "a2_median":np.median(anchor_2_counts),
        "a2_ci":    ci_2,
        "a2_p":     p_anchor_2,
        # Month of 80th event
        "m80_mean":   finite_80.mean() if len(finite_80) > 0 else np.nan,
        "m80_median": np.median(finite_80) if len(finite_80) > 0 else np.nan,
        "m80_ci":     np.percentile(finite_80, [2.5, 97.5]) if len(finite_80) > 0 else [np.nan, np.nan],
        "pct_reached_80": len(finite_80) / len(month_80_vals) * 100,
    }


# ─────────────────────────────────────────────
# Cure Fraction Sweep
# ─────────────────────────────────────────────

def sweep_cure_fractions(
    p: ModelParams,
    timeline: TrialTimeline,
    cure_fractions: np.ndarray,
    n_sims: int = 5_000,
    seed: int = 42
) -> dict:
    """
    Sweep across a range of cure fractions and find which values
    are consistent with both SEC anchor observations.
    """
    a1_means, a1_cis_lo, a1_cis_hi = [], [], []
    a2_means, a2_cis_lo, a2_cis_hi = [], [], []

    for cf in cure_fractions:
        print(f"  π={cf*100:.0f}%", end=" ")
        mc = run_monte_carlo(p, timeline, cf, n_sims=n_sims, seed=seed)
        a1_means.append(mc["a1_mean"])
        a1_cis_lo.append(mc["a1_ci"][0])
        a1_cis_hi.append(mc["a1_ci"][1])
        a2_means.append(mc["a2_mean"])
        a2_cis_lo.append(mc["a2_ci"][0])
        a2_cis_hi.append(mc["a2_ci"][1])

    return {
        "cure_fractions": cure_fractions,
        "a1_means":   np.array(a1_means),
        "a1_cis_lo":  np.array(a1_cis_lo),
        "a1_cis_hi":  np.array(a1_cis_hi),
        "a2_means":   np.array(a2_means),
        "a2_cis_lo":  np.array(a2_cis_lo),
        "a2_cis_hi":  np.array(a2_cis_hi),
    }


# ─────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────

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


def plot_simulation_results(
    mc: dict,
    sweep: dict,
    p: ModelParams,
    timeline: TrialTimeline,
    cure_fraction: float,
    single_trial: dict
):
    fig = plt.figure(figsize=(18, 13), facecolor=DARK_BG)
    fig.suptitle(
        "GPS · REGAL  |  Monte Carlo Enrollment Simulation",
        color=TEXT_HI, fontsize=16, fontfamily="monospace", y=0.98, x=0.02, ha="left"
    )
    fig.text(0.02, 0.955,
             f"n=10,000 simulations  ·  π={cure_fraction*100:.1f}%  ·  "
             f"Enrollment Jan 2020 – Jun 2022  ·  Anchors: 60@mo46, 72@mo58",
             color=TEXT_DIM, fontsize=9, fontfamily="monospace")

    gs = gridspec.GridSpec(3, 3, figure=fig,
                           hspace=0.45, wspace=0.35,
                           left=0.06, right=0.97, top=0.93, bottom=0.06)

    # ── Stat banner ──────────────────────────────────────────────────
    obs_in_ci_1 = mc["a1_ci"][0] <= timeline.anchor_1_events <= mc["a1_ci"][1]
    obs_in_ci_2 = mc["a2_ci"][0] <= timeline.anchor_2_events <= mc["a2_ci"][1]

    stats = [
        ("Cure Fraction π",    f"{cure_fraction*100:.1f}%",              BLUE),
        ("Anchor 1 (mo46)",    f"obs={timeline.anchor_1_events} | sim={mc['a1_mean']:.1f}", GREEN if obs_in_ci_1 else RED),
        ("Anchor 2 (mo58)",    f"obs={timeline.anchor_2_events} | sim={mc['a2_mean']:.1f}", GREEN if obs_in_ci_2 else RED),
        ("95% CI mo46",        f"[{mc['a1_ci'][0]:.0f}, {mc['a1_ci'][1]:.0f}]",            GREEN if obs_in_ci_1 else RED),
        ("95% CI mo58",        f"[{mc['a2_ci'][0]:.0f}, {mc['a2_ci'][1]:.0f}]",            GREEN if obs_in_ci_2 else RED),
        ("80th event (median)",f"mo {mc['m80_median']:.1f}",              AMBER),
    ]
    for idx, (lbl, val, col) in enumerate(stats):
        x = 0.02 + idx * 0.163
        fig.text(x, 0.66, lbl, color=TEXT_DIM, fontsize=7,
                 fontfamily="monospace", transform=fig.transFigure)
        fig.text(x, 0.64, val, color=col, fontsize=11,
                 fontfamily="monospace", fontweight="bold",
                 transform=fig.transFigure)

    # ── Panel 1: Enrollment timeline (single sim) ─────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    gps_mask = single_trial["arms"] == "GPS"
    ax1.hist(single_trial["enrollment_times"][gps_mask],
             bins=20, color=BLUE, alpha=0.7, label="GPS")
    ax1.hist(single_trial["enrollment_times"][~gps_mask],
             bins=20, color=ORANGE, alpha=0.7, label="BAT")
    ax1.axvline(timeline.enrollment_end, color=TEXT_DIM, ls="--", lw=1,
                label="Accrual close")
    ax1.set_xlabel("Calendar Month from Trial Start")
    ax1.set_ylabel("Patients Enrolled")
    ax1.legend(fontsize=7, facecolor=DARK_BG, edgecolor=GRID_COL,
               labelcolor=TEXT_MID)
    style_ax(ax1, "ENROLLMENT DISTRIBUTION (1 sim)")

    # ── Panel 2: Cumulative events over calendar time (single sim) ────
    ax2 = fig.add_subplot(gs[0, 1])

    cal_times = single_trial["calendar_event_times"]
    finite_mask = np.isfinite(cal_times)
    sorted_cal = np.sort(cal_times[finite_mask & single_trial["had_event"]])

    # GPS and BAT separately
    gps_events = np.sort(cal_times[gps_mask & single_trial["had_event"] & finite_mask])
    bat_events = np.sort(cal_times[~gps_mask & single_trial["had_event"] & finite_mask])

    if len(gps_events) > 0:
        ax2.step(gps_events, np.arange(1, len(gps_events)+1),
                 color=BLUE, lw=2, where="post", label="GPS events")
    if len(bat_events) > 0:
        ax2.step(bat_events, np.arange(1, len(bat_events)+1),
                 color=ORANGE, lw=2, where="post", label="BAT events")

    # Total pooled
    if len(sorted_cal) > 0:
        ax2.step(sorted_cal, np.arange(1, len(sorted_cal)+1),
                 color=PURPLE, lw=1.5, where="post", ls="--", label="Pooled")

    # Anchor lines
    ax2.axvline(timeline.anchor_1_month, color=GREEN, ls=":", lw=1.5,
                label=f"mo{timeline.anchor_1_month:.0f}: obs={timeline.anchor_1_events}")
    ax2.axvline(timeline.anchor_2_month, color=AMBER, ls=":", lw=1.5,
                label=f"mo{timeline.anchor_2_month:.0f}: obs={timeline.anchor_2_events}")
    ax2.axhline(timeline.anchor_1_events, color=GREEN, ls=":", lw=0.8, alpha=0.5)
    ax2.axhline(timeline.anchor_2_events, color=AMBER, ls=":", lw=0.8, alpha=0.5)
    ax2.axhline(timeline.final_analysis_events, color=RED, ls="--", lw=1,
                label=f"Final analysis ({timeline.final_analysis_events} events)")

    ax2.set_xlabel("Calendar Month from Trial Start")
    ax2.set_ylabel("Cumulative Events")
    ax2.set_xlim(0, timeline.calendar_cutoff)
    ax2.legend(fontsize=6, facecolor=DARK_BG, edgecolor=GRID_COL,
               labelcolor=TEXT_MID, loc="upper left")
    style_ax(ax2, "CUMULATIVE EVENTS OVER TIME (1 sim)")

    # ── Panel 3: Distribution of events at anchor 1 ───────────────────
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.hist(mc["anchor_1_counts"], bins=30, color=BLUE, alpha=0.8,
             edgecolor=DARK_BG, linewidth=0.3)
    ax3.axvline(timeline.anchor_1_events, color=GREEN, lw=2.5,
                label=f"Observed: {timeline.anchor_1_events}")
    ax3.axvline(mc["a1_mean"], color=BLUE, lw=1.5, ls="--",
                label=f"Sim mean: {mc['a1_mean']:.1f}")
    ax3.axvline(mc["a1_ci"][0], color=TEXT_DIM, lw=1, ls=":",
                label=f"95% CI: [{mc['a1_ci'][0]:.0f}, {mc['a1_ci'][1]:.0f}]")
    ax3.axvline(mc["a1_ci"][1], color=TEXT_DIM, lw=1, ls=":")
    ax3.set_xlabel("Events at Calendar Month 46")
    ax3.set_ylabel("Simulation Count")
    ax3.legend(fontsize=7, facecolor=DARK_BG, edgecolor=GRID_COL,
               labelcolor=TEXT_MID)
    consistency_1 = "✓ CONSISTENT" if obs_in_ci_1 else "✗ INCONSISTENT"
    ax3.set_title(f"ANCHOR 1 VALIDATION  ·  {consistency_1}",
                  color=GREEN if obs_in_ci_1 else RED,
                  fontsize=9, fontfamily="monospace", pad=10, loc="left")
    for spine in ax3.spines.values():
        spine.set_edgecolor(GRID_COL)
    ax3.set_facecolor(PANEL_BG)
    ax3.tick_params(colors=TEXT_DIM, labelsize=8)
    ax3.grid(True, color=GRID_COL, linewidth=0.5, linestyle="--", alpha=0.7)

    # ── Panel 4: Distribution of events at anchor 2 ───────────────────
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.hist(mc["anchor_2_counts"], bins=30, color=ORANGE, alpha=0.8,
             edgecolor=DARK_BG, linewidth=0.3)
    ax4.axvline(timeline.anchor_2_events, color=GREEN, lw=2.5,
                label=f"Observed: {timeline.anchor_2_events}")
    ax4.axvline(mc["a2_mean"], color=ORANGE, lw=1.5, ls="--",
                label=f"Sim mean: {mc['a2_mean']:.1f}")
    ax4.axvline(mc["a2_ci"][0], color=TEXT_DIM, lw=1, ls=":")
    ax4.axvline(mc["a2_ci"][1], color=TEXT_DIM, lw=1, ls=":",
                label=f"95% CI: [{mc['a2_ci'][0]:.0f}, {mc['a2_ci'][1]:.0f}]")
    ax4.set_xlabel("Events at Calendar Month 58")
    ax4.set_ylabel("Simulation Count")
    ax4.legend(fontsize=7, facecolor=DARK_BG, edgecolor=GRID_COL,
               labelcolor=TEXT_MID)
    consistency_2 = "✓ CONSISTENT" if obs_in_ci_2 else "✗ INCONSISTENT"
    ax4.set_title(f"ANCHOR 2 VALIDATION  ·  {consistency_2}",
                  color=GREEN if obs_in_ci_2 else RED,
                  fontsize=9, fontfamily="monospace", pad=10, loc="left")
    for spine in ax4.spines.values():
        spine.set_edgecolor(GRID_COL)
    ax4.set_facecolor(PANEL_BG)
    ax4.tick_params(colors=TEXT_DIM, labelsize=8)
    ax4.grid(True, color=GRID_COL, linewidth=0.5, linestyle="--", alpha=0.7)

    # ── Panel 5: Distribution of month of 80th event ──────────────────
    ax5 = fig.add_subplot(gs[1, 1])
    finite_80 = mc["finite_80"]
    if len(finite_80) > 0:
        ax5.hist(finite_80, bins=30, color=PURPLE, alpha=0.8,
                 edgecolor=DARK_BG, linewidth=0.3)
        ax5.axvline(mc["m80_median"], color=AMBER, lw=2.5,
                    label=f"Median: mo {mc['m80_median']:.1f}")
        ax5.axvline(mc["m80_ci"][0], color=TEXT_DIM, lw=1, ls=":")
        ax5.axvline(mc["m80_ci"][1], color=TEXT_DIM, lw=1, ls=":",
                    label=f"95% CI: mo [{mc['m80_ci'][0]:.1f}, {mc['m80_ci'][1]:.1f}]")
        # Current date marker (mo ~72 = Apr 2026)
        ax5.axvline(72, color=GREEN, lw=2, ls="--", label="Now (~mo 72)")
        ax5.legend(fontsize=7, facecolor=DARK_BG, edgecolor=GRID_COL,
                   labelcolor=TEXT_MID)
    ax5.set_xlabel("Calendar Month of 80th Event")
    ax5.set_ylabel("Simulation Count")
    style_ax(ax5, f"FINAL ANALYSIS TIMING  ·  {mc['pct_reached_80']:.1f}% sims reach 80 events")

    # ── Panel 6: Cure fraction sweep ──────────────────────────────────
    ax6 = fig.add_subplot(gs[1, 2])
    cf_pct = sweep["cure_fractions"] * 100

    ax6.fill_between(cf_pct, sweep["a1_cis_lo"], sweep["a1_cis_hi"],
                     alpha=0.25, color=BLUE, label="95% CI (mo46)")
    ax6.plot(cf_pct, sweep["a1_means"], color=BLUE, lw=2, label="Sim mean (mo46)")

    ax6.fill_between(cf_pct, sweep["a2_cis_lo"], sweep["a2_cis_hi"],
                     alpha=0.25, color=ORANGE, label="95% CI (mo58)")
    ax6.plot(cf_pct, sweep["a2_means"], color=ORANGE, lw=2, label="Sim mean (mo58)")

    # Observed anchor lines
    ax6.axhline(timeline.anchor_1_events, color=BLUE, ls="--", lw=1.5,
                label=f"Observed mo46={timeline.anchor_1_events}")
    ax6.axhline(timeline.anchor_2_events, color=ORANGE, ls="--", lw=1.5,
                label=f"Observed mo58={timeline.anchor_2_events}")

    # Mark base cure fraction
    ax6.axvline(cure_fraction * 100, color=GREEN, ls=":", lw=1.5,
                label=f"Base π={cure_fraction*100:.1f}%")

    ax6.set_xlabel("Cure Fraction π (%)")
    ax6.set_ylabel("Expected Events")
    ax6.legend(fontsize=6, facecolor=DARK_BG, edgecolor=GRID_COL,
               labelcolor=TEXT_MID, loc="upper right")
    style_ax(ax6, "CURE FRACTION SWEEP  ·  Which π fits the anchors?")

    # ── Panel 7: Follow-up distribution (single sim) ──────────────────
    ax7 = fig.add_subplot(gs[2, :2])
    gps_fu = single_trial["follow_up"][gps_mask]
    bat_fu = single_trial["follow_up"][~gps_mask]

    bins = np.linspace(0, timeline.calendar_cutoff, 25)
    ax7.hist(bat_fu, bins=bins, color=ORANGE, alpha=0.7, label="BAT follow-up")
    ax7.hist(gps_fu, bins=bins, color=BLUE, alpha=0.7, label="GPS follow-up")

    # Overlay: who had event vs censored
    gps_event_fu = single_trial["follow_up"][gps_mask & single_trial["had_event"]]
    bat_event_fu = single_trial["follow_up"][~gps_mask & single_trial["had_event"]]
    ax7.hist(gps_event_fu, bins=bins, color=BLUE, alpha=1.0,
             histtype="step", lw=2, label="GPS events (uncensored)")
    ax7.hist(bat_event_fu, bins=bins, color=ORANGE, alpha=1.0,
             histtype="step", lw=2, label="BAT events (uncensored)")

    ax7.set_xlabel("Follow-up Time from Enrollment (months)")
    ax7.set_ylabel("Patients")
    ax7.legend(fontsize=7, facecolor=DARK_BG, edgecolor=GRID_COL,
               labelcolor=TEXT_MID)
    style_ax(ax7, "FOLLOW-UP DISTRIBUTION  ·  Filled=all patients, Outline=had event (1 sim)")

    # ── Panel 8: Summary text ─────────────────────────────────────────
    ax8 = fig.add_subplot(gs[2, 2])
    ax8.set_facecolor(PANEL_BG)
    ax8.axis("off")
    for spine in ax8.spines.values():
        spine.set_edgecolor(GRID_COL)

    lines = [
        ("VALIDATION SUMMARY", TEXT_MID, 10, True),
        ("", TEXT_DIM, 8, False),
        (f"Anchor 1 (mo46): obs={timeline.anchor_1_events}", TEXT_DIM, 8, False),
        (f"  sim mean: {mc['a1_mean']:.1f}", BLUE, 9, False),
        (f"  95% CI:   [{mc['a1_ci'][0]:.0f}, {mc['a1_ci'][1]:.0f}]", BLUE, 9, False),
        (f"  → {consistency_1}", GREEN if obs_in_ci_1 else RED, 9, True),
        ("", TEXT_DIM, 8, False),
        (f"Anchor 2 (mo58): obs={timeline.anchor_2_events}", TEXT_DIM, 8, False),
        (f"  sim mean: {mc['a2_mean']:.1f}", ORANGE, 9, False),
        (f"  95% CI:   [{mc['a2_ci'][0]:.0f}, {mc['a2_ci'][1]:.0f}]", ORANGE, 9, False),
        (f"  → {consistency_2}", GREEN if obs_in_ci_2 else RED, 9, True),
        ("", TEXT_DIM, 8, False),
        (f"80th event:", TEXT_DIM, 8, False),
        (f"  median: mo {mc['m80_median']:.1f}", PURPLE, 9, False),
        (f"  95% CI: [{mc['m80_ci'][0]:.1f}, {mc['m80_ci'][1]:.1f}]", PURPLE, 9, False),
        (f"  % sims reaching 80: {mc['pct_reached_80']:.1f}%", PURPLE, 9, False),
        ("", TEXT_DIM, 8, False),
        ("Enrollment: uniform Jan2020–Jun2022", TEXT_DIM, 7, False),
        (f"n_sims: 10,000  ·  seed: 42", TEXT_DIM, 7, False),
    ]
    y = 0.95
    for text, color, size, bold in lines:
        ax8.text(0.05, y, text, color=color, fontsize=size,
                 fontfamily="monospace",
                 fontweight="bold" if bold else "normal",
                 transform=ax8.transAxes)
        y -= 0.055

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "gps_enrollment_sim_output.png"
    )
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=DARK_BG)
    print(f"  Saved: {out_path}")
    plt.show()


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

if __name__ == "__main__":
    p        = ModelParams()
    timeline = TrialTimeline()

    # ── Step 1: Sweep first to find implied π ──────────────────────
    print("\n  [1/2] Sweeping cure fractions (5,000 sims each)...")
    cf_range = np.linspace(0.05, 0.45, 9)
    sweep = sweep_cure_fractions(p, timeline, cf_range, n_sims=5_000, seed=42)
    print()

    cf_pct = sweep["cure_fractions"] * 100

    # Interpolate — note: means are DECREASING as π increases
    # so we need to sort for interp1d to work correctly
    a1_means = sweep["a1_means"]
    a2_means = sweep["a2_means"]

    # Sort by ascending mean (interp1d needs monotonic x)
    sort_idx_1 = np.argsort(a1_means)
    sort_idx_2 = np.argsort(a2_means)

    f_anchor1 = interp1d(
        a1_means[sort_idx_1], cf_pct[sort_idx_1],
        kind="linear", fill_value="extrapolate"
    )
    f_anchor2 = interp1d(
        a2_means[sort_idx_2], cf_pct[sort_idx_2],
        kind="linear", fill_value="extrapolate"
    )

    implied_pi_1 = float(f_anchor1(timeline.anchor_1_events))
    implied_pi_2 = float(f_anchor2(timeline.anchor_2_events))
    needed_cure_fraction = (implied_pi_1 + implied_pi_2) / 2 / 100  # convert % back to fraction

    print(f"  Implied π from anchor 1 (mo46): {implied_pi_1:.1f}%")
    print(f"  Implied π from anchor 2 (mo58): {implied_pi_2:.1f}%")
    print(f"  Average implied π:              {needed_cure_fraction*100:.1f}%")

    # ── Step 2: Everything else uses needed_cure_fraction ──────────
    print("=" * 60)
    print("  GPS · REGAL · Monte Carlo Enrollment Simulation")
    print("=" * 60)
    print(f"  Trial start:          Jan 2020 (month 0)")
    print(f"  Enrollment close:     ~Jun 2022 (month {timeline.enrollment_end:.0f})")
    print(f"  Calendar cutoff:      ~Jan 2026 (month {timeline.calendar_cutoff:.0f})")
    print(f"  Anchor 1:             {timeline.anchor_1_events} events at month {timeline.anchor_1_month:.0f}")
    print(f"  Anchor 2:             {timeline.anchor_2_events} events at month {timeline.anchor_2_month:.0f}")
    print(f"  Implied cure fraction π: {needed_cure_fraction*100:.1f}%")
    print(f"  BAT median OS:        {p.bat_median_os}mo")
    print(f"  GPS susceptible mOS:  {p.gps_median_os}mo")
    print("-" * 60)

    # Single trial visualization
    rng_single = np.random.default_rng(99)
    single = simulate_trial(p, timeline, needed_cure_fraction, rng_single)
    print(f"  Single sim — events at mo46: {single['n_at_anchor_1']}, "
          f"mo58: {single['n_at_anchor_2']}, "
          f"80th event: mo {single['month_of_80th_event']:.1f}")

    # ── Step 3: Main Monte Carlo with implied π ─────────────────────
    print("\n  [2/2] Running main Monte Carlo (10,000 sims)...")
    mc = run_monte_carlo(p, timeline, needed_cure_fraction, n_sims=10_000, seed=42)

    print(f"\n  Anchor 1 (mo46):")
    print(f"    Observed:  {timeline.anchor_1_events}")
    print(f"    Sim mean:  {mc['a1_mean']:.1f}")
    print(f"    95% CI:    [{mc['a1_ci'][0]:.0f}, {mc['a1_ci'][1]:.0f}]")
    obs_in_1 = mc["a1_ci"][0] <= timeline.anchor_1_events <= mc["a1_ci"][1]
    print(f"    Consistent: {'YES ✓' if obs_in_1 else 'NO ✗'}")

    print(f"\n  Anchor 2 (mo58):")
    print(f"    Observed:  {timeline.anchor_2_events}")
    print(f"    Sim mean:  {mc['a2_mean']:.1f}")
    print(f"    95% CI:    [{mc['a2_ci'][0]:.0f}, {mc['a2_ci'][1]:.0f}]")
    obs_in_2 = mc["a2_ci"][0] <= timeline.anchor_2_events <= mc["a2_ci"][1]
    print(f"    Consistent: {'YES ✓' if obs_in_2 else 'NO ✗'}")

    print(f"\n  80th event (final analysis):")
    print(f"    Median calendar month: {mc['m80_median']:.1f}")
    print(f"    95% CI: [{mc['m80_ci'][0]:.1f}, {mc['m80_ci'][1]:.1f}]")
    print(f"    % sims reaching 80 events: {mc['pct_reached_80']:.1f}%")

    # Plot — pass needed_cure_fraction everywhere
    plot_simulation_results(mc, sweep, p, timeline, needed_cure_fraction, single)
    print("=" * 60)
