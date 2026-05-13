---
type: summary
title: "$SLS Deepest Due Diligence for REGAL Trial (From a Deep Value Investor)"
source_type: reddit
author: Confident-Web-7118
subreddit: r/pennystocks
published: 2026-02-15
date_ingested: 2026-04-11
original_file: "raw/reddit/$SLS (Deepest Due Diligence for REGAL Trial) (From a Deep Value Investor).md"
tags: [regal, gps, cure-fraction, event-velocity, reverse-engineering, idmc, bat-arm, retail-dd, conflict-of-interest]
conflict_of_interest: "Author holds 805,000 shares of SELS and is continuously accumulating"
---

## Summary

A detailed retail investor DD post from Feb 2026 by Confident-Web-7118 — the same author whose earlier (Apr 2025) r/pennystocks post was critiqued in the prior Claude session. This post is substantially more rigorous: it presents a reverse-engineered cure fraction model anchored to verified public data (72 events as of Dec 26, 2025), a velocity proof, enrollment S-curve stress testing, and an explicit devil's advocate section. The author discloses holding 805K shares.

The post adds critical new facts not previously in the wiki: the trial N is confirmed at 126, the IDMC has issued "continue without modification" at **both** interim reviews, the enrollment period is now defined (Feb 2021 – Apr 2024), and the velocity proof provides the strongest publicly available quantitative evidence for a cure fraction plateau.

**Methodological quality: Medium.** Substantially better than the earlier post. The core arguments (velocity deceleration, BAT floor, impossible-to-fail sensitivity analysis) are genuine and defensible. The specific numbers (42–48% cure fraction, "99.9999% probability") remain retail-grade inference from aggregate data, not patient-level analysis.

---

## Key New Facts (not previously in wiki)

### Trial Demographics
- **N = 126** (confirmed; not 140 as originally planned; ~63 per arm)
- **Enrollment period:** February 2021 (first patient) → April 2024 (126th patient)
- **Protocol amendment: November 2022** — likely accelerated late enrollment (S-curve midpoint ~month 25, consistent with this)
- **Event maturity:** 72/126 = 57.1% — past the pooled median OS (a "hard historical fact")

### Event Count & Timing (PUBLIC DATA)
- **60 events by December 2024** (consistent with mo46 SEC anchor)
- **72 events by December 26, 2025** (12.5 months later)
- **54 patients still alive** at Dec 26, 2025
- **66 patients at risk** at the Dec 2024 → Dec 2025 window start

### Velocity Proof (Strongest Evidence)
- 12 deaths over 12.5 months from 66 patients at risk
- Hazard rate: 12 ÷ (66 × 12.5) = **0.0145 deaths/person-month**
- Annualized mortality: **16%**
- Implied median survival for the surviving population: **~48 months**

Comparison table (expected deaths from 66 patients over 12.5 months at various mOS):

| mOS assumption | Expected deaths |
|---|---|
| 10 months | 38.3 |
| 14.5 months | 29.7 |
| 20 months | 23.2 |
| 30 months | 16.6 |
| 50 months | 10.5 |
| **OBSERVED** | **12** |

The observed rate matches a surviving population with implied mOS ~48 months. This is the strongest public evidence of a cure fraction plateau — the survival curve has effectively flatlined.

### IDMC: Two "Continue Without Modification" Rulings
- IDMC recommended "continue without modification" at **both interim reviews**
- This is an update from the wiki's prior knowledge (one known interim at 60 events)
- Absence of futility stop at both reviews is a meaningful positive signal
- Note: O'Brien-Fleming boundaries at interim are conservative; continuation is common but a futility stop would have been definitively negative

### Arm Survival Breakdown (Model Output, Not Official Data)
Author's model predicts the 54 surviving patients break down as:

| | BAT Arm | GPS Arm |
|---|---|---|
| Total enrolled | 63 | 63 |
| Projected dead | ~57 | ~18 |
| **Projected alive** | **~6 (10%)** | **~45 (71%)** |
| "Cured" (GPS plateau) | — | ~26–30 |

*Caveat: these are model outputs, not official data.*

### Phase 2 CR2 (Brayer/Moffitt — NEW SOURCE)
Author references a Phase 2 GPS study specifically in **CR2 patients** (same population as REGAL):
- **GPS mOS: 21.0 months**
- **Control mOS: 5.4 months**
- **No cure fraction plateau observed** — attributed to fixed dosing (6 shots then stop)
- This is distinct from Maslak 2018 (which appears to be primarily CR1)

See [[entities/gps]] for the dosing-change hypothesis.

### Dosing Change Hypothesis (Key Mechanistic Argument)
Author argues this is the critical REGAL design improvement:

| Feature | Phase 2 CR2 | Phase 3 REGAL |
|---|---|---|
| Dosing | 6 shots, then stop | Monthly boosters indefinitely |
| Duration | Fixed schedule | Until relapse |
| Observed mOS (GPS) | 21.0 months | Modeled >60+ months |
| Cure fraction | None observed | 42–48% modeled |

Hypothesis: continuous boosting in REGAL converts "delayed death" (Phase 2 CR2) into "long-term immune surveillance" — reproducing the CR1 ghost plateau (47%) in a CR2 population.

### Expected Readout Timeline
- **80th event (final trigger):** Q2–Q3 2026
- The trial may never hit 80 events organically (asymptotic max ~93 if plateau is real)
- SELLAS may trigger final analysis on a calendar date rather than waiting for the 80th event

---

## Model Claims (With Caveats)

### Reverse-Engineered Cure Fraction: 42–48%
- **Method:** Constrained by the 72-event count at month 58 + enrollment S-curve
- **Base case (BAT=10mo):** 42% cure fraction, uncured GPS mOS 34–39mo, GPS theoretical mOS 97–183 months
- **Unconstrained grid search:** BAT floats to 14.5mo, cure fraction still 64%
- **Conservative (back-loaded enrollment):** BAT drops to 12.5–13mo, cure fraction stays at 64%

**Key discrepancy with our model:** Our biology-derived π = 29.5%; author's reverse-engineered π = 42–48%. The gap is significant and likely reflects:
1. Author assumes ALL surviving patients are on one curve; we separately modeled HLA subgroups
2. The Phase 2 dosing change (continuous vs. fixed) may genuinely have produced a higher cure fraction than Phase 2 biology predicts
3. Upward bias in retail investor reverse-engineering from summary data
4. Our model may underestimate cure fraction by not accounting for non-A*02:01 HLA coverage (see [[sources/doubrovina-2012-wt1-hla-mapping]])

### Expected Topline HR: 0.35–0.50
- Conditional HR (responders only at BAT=10mo): **~0.13** — mathematically correct for the cure fraction subpopulation
- Adjusted for early non-responder deaths dragging the Cox average: **0.31–0.49**
- REGAL threshold: HR < 0.636

Sensitivity table:

| BAT mOS | Conditional HR (responders) | P(success) |
|---|---|---|
| 8 mo | 0.10 | 100% |
| 10 mo | 0.13 | 100% |
| 12 mo | 0.16 | 100% |
| 14 mo | 0.22 | 100% |
| 16 mo | 0.31 | 100% |
| 18 mo | 0.45 | ~99% |
| 20 mo | 0.61 | ~95% |

**Trial fails only if BAT mOS > 23 months** — author states this has never been observed in AML CR2 (not eligible for transplant) history.

---

## Critical Assessment

### What This Post Gets Right
- **Velocity proof is genuinely compelling.** 12 deaths/12.5 months from 66 at risk is a real, verifiable calculation anchored to public data. The implied mOS of ~48 months for the surviving cohort cannot be explained by a standard Weibull model without a cure fraction.
- **BAT failure floor analysis is solid.** The 23-month BAT mOS required for trial failure exceeds any documented historical benchmark by 5+ months. Combined with allo-SCT exclusion, this is a genuine margin of safety.
- **Enrollment S-curve stress testing is methodologically appropriate.**
- **Disclosure of conflict of interest** (805K shares, continuous accumulation).

### What This Post Gets Wrong or Overstates
- **"99.9999% probability"**: Still present in the Stocktwits excerpts. Indefensible from aggregate summary data. Should be read as "very confident", not as a statistical claim.
- **42–48% cure fraction is reverse-engineered from two anchor points**: The underlying model has many degrees of freedom. The "must be cure fraction" argument is directionally correct but the specific number could vary widely.
- **Survivor arm breakdown (~45 GPS alive) is a model output**: Not verified official data. Could be directionally right.
- **Phase 2 CR2 Brayer/Moffitt data cited without a source**: This is a new and important data point (21.0 vs 5.4 months) but the original paper is not cited. Needs verification.

---

## Relevance to Outcome Prediction

This post materially strengthens the bull case by providing the most current publicly available event velocity data (Dec 2025). The velocity proof alone — 12 deaths from 66 at risk over 12.5 months — is the strongest evidence yet of a cure fraction plateau in real-time REGAL data.

The new IDMC "continue without modification" (both reviews) rules out a futility stop having occurred, which was a key residual risk.

The Phase 2 CR2 dosing hypothesis (continuous boosters vs. fixed 6 shots) is the single most important mechanistic argument in this post — it explains why REGAL might show a cure fraction that the Phase 2 CR2 did not.

## Open Questions Raised

- What is the source for the Phase 2 CR2 Brayer/Moffitt data (GPS mOS 21.0 vs control 5.4 months)?
- Was there truly a second IDMC review beyond the 60-event interim, or is the author conflating two reviews of the same interim?
- What triggered the November 2022 protocol amendment, and did it change anything beyond enrollment practices?
- Has SELLAS communicated any calendar-date trigger plan for final analysis?
