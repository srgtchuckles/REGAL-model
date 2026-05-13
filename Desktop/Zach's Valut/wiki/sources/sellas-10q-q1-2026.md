---
type: source
title: "SELLAS Life Sciences Group — Form 10-Q, Q1 2026 (Quarter Ended March 31, 2026)"
source_type: sec-filing
date_ingested: 2026-05-12
original_file: wiki/sources/sellas-10q-q1-2026.pdf
filed: 2026-05-12
period: Q1 2026 (ended March 31, 2026)
tags: [sellas, regal, gps, event-count, bla-prep, financials, sec-filing, sls009, 3dmed-arbitration]
---

## Summary

SELLAS's quarterly SEC filing for Q1 2026, filed May 12, 2026. The document is primarily a financial filing, but contains the most current publicly available event count for the REGAL trial and material signals about management's internal expectations.

**The single most important fact:** As of May 11, 2026 (the day before filing), the pooled number of REGAL events was **78 of 80** required to trigger final analysis. The trial is 2 deaths away from readout.

The filing also reveals an aggressive BLA preparation spend in Q1 2026 — $1.3M increase in manufacturing costs and $0.4M increase in regulatory consulting, both explicitly attributed to "preparing for a potential BLA filing following final analysis of the REGAL study." This spending pattern is a meaningful secondary signal of management confidence, though it carries corporate incentive bias.

---

## Key Claims

### REGAL Trial — Event Count (CRITICAL)
- **78 events as of May 11, 2026** (per CRO report to company)
- Final analysis triggered at **80 events** (deaths)
- Company remains fully blinded to all efficacy and survival data
- Filing explicitly states: "this one-time update on the aggregate number of events does not impact future statistical analyses" (i.e., no statistical penalty for this disclosure)
- "We will announce the 80th event when it occurs"

### REGAL Trial — Timeline & Status
- Enrollment: 125–140 patients enrolled at ~95 clinical sites; March 2024 completion confirmed
- Interim analysis: Triggered December 2024 (60 events). IDMC completed in **January 2025** and recommended "continue without modifications"
- No second interim analysis mentioned; "the next and final analysis will be conducted once 80 events are reached"
- The filing does **not** announce the 80th event, so the final analysis had not been triggered as of May 12, 2026

### BLA Preparation Spend (Management Signal)
Q1 2026 R&D expenses breakdown for GPS:
- External clinical trial expenses (GPS): $987K (vs $741K in Q1 2025, +$246K)
- Manufacturing and clinical drug supply: $1,335K total R&D line — **+$1,257K YoY**, explicitly "to prepare for a potential BLA filing"
- Clinical and regulatory consulting: +$399K YoY, explicitly for BLA prep

This level of BLA prep spend before data readout is consistent with management anticipating a positive result, though all pre-commercialization spend is speculative.

### Financial Position
- Cash and equivalents: **$107.1M** as of March 31, 2026
- Post-Q1 additional proceeds: $7.5M from 4.7M warrant exercises at avg $1.60/share
- **Total effective cash: ~$114.6M** as of filing date
- Accumulated deficit: $283.4M
- 2026 ATM facility: up to $150M with TD Cowen (established March 2026, not yet used)
- Q1 2026 net loss: $8.4M ($9.2M opex, $0.8M interest income)
- Cash runway: "at least the next twelve months from the date of issuance" (i.e., no going concern risk)
- Shares outstanding: **186,032,574** as of May 11, 2026
- Stock price reference: $4.23/share on March 31, 2026 (NASDAQ close)

### SLS009 (Tambiciclib) — Secondary Program Status
- Phase 2 trial in r/r AML met all primary endpoints (July 2025 announcement referenced)
- ORR 33% overall; 40% at 30mg BIW dose; 44% in AML MR subgroup
- Median OS 8.9 months in AML MR patients; 8.8 months in venetoclax-refractory patients at 30mg BIW
- FDA recommended advancing to first-line study
- New randomized 80-patient Phase 2 trial began enrollment Q1 2026 (newly diagnosed + early venetoclax resistance cohorts)
- IMPACT-AML European partnership announced January 2026

### 3D Medicines Arbitration
- Binding arbitration commenced December 2023 at HKIAC
- Evidence hearing held **January 2026** (HKIAC)
- Final determination pending as of March 31, 2026
- Dispute: milestone payments ($13M claimed) + failure to use commercially reasonable best efforts in China
- $191.5M in potential future milestones still outstanding under 3D Medicines Agreement

---

## Relevance to REGAL Outcome Prediction

### Event Count Is the Most Important Update in the Wiki

The 78/80 event data point, dated May 11, 2026, is the most current and precise public anchor available. Combined with the prior anchor from [[sources/confident-web-7118-regal-dd-feb2026]] (72 events as of Dec 26, 2025), we can now reconstruct the trailing event velocity:

| Period | Events | Duration | Rate |
|---|---|---|---|
| Dec 2024 → Dec 26, 2025 | 60 → 72 (+12) | ~12 months | 1.0 events/mo |
| Dec 26, 2025 → May 11, 2026 | 72 → 78 (+6) | ~4.5 months | 1.33 events/mo |

The rate from Dec 2025 → May 2026 is **1.33 events/month** from a surviving cohort of ~48 patients (54 alive at Dec 2025 minus ~6 subsequent deaths). At this rate, the remaining 2 events would arrive in **approximately 1.5 months** from May 11, placing the expected final analysis trigger around **late June – early July 2026**.

However, this is subject to high variance (only 2 events needed, from ~48 survivors). The actual trigger could be sooner or later.

### Velocity Interpretation
The deceleration from the earlier period (pre-Dec 2024) to now is consistent with a cure fraction plateau — the faster events early in the trial represent the non-immune-responder population dying, while the survivors represent a slower-dying or non-dying cohort. A rate of 1.33 events/month from ~48 patients implies an annualized mortality rate of ~33%/year for the surviving cohort, which is higher than the prior Claude session's most optimistic scenario but still far below what a simple Weibull model would predict without a cure fraction.

### BLA Prep as a Management Signal
The spending ramp-up ($1.3M manufacturing, $0.4M consulting) ahead of data is consistent with, but does not prove, management confidence. Companies do pre-position for FDA submissions regardless of outcome to minimize commercialization lag if positive. However, the specificity and scale of Q1 2026 spend relative to Q1 2025 is notable. This is a weak positive signal, not strong evidence.

### Cash Position
$114.6M in effective cash is sufficient to fund operations through a BLA submission and FDA review if the trial succeeds. No emergency dilution risk. This is a background supportive condition, not predictive.

---

## Contradictions / Tensions

- No tension with other wiki sources. This filing confirms the IDMC January 2025 "continue without modification" signal already known from other sources.
- The event velocity from Dec 2025 → May 2026 (1.33/month) is slightly higher than the extreme optimistic plateau scenario would predict, but not inconsistent with a mixture of cured and non-cured survivors. Needs no model revision.

---

## Open Questions

- When exactly will the 80th event be announced? (Company has committed to immediate disclosure)
- Will the final analysis produce OS data only, or also EFS/PFS secondary endpoints?
- Is there any calendar-date trigger provision if the 80th event takes longer than expected? (Not mentioned in this filing)
- 3DMed arbitration outcome — if the $13M milestone is awarded, does it affect SELLAS cash position materially? ($13M vs $114M cash = 11%, notable but not critical)
