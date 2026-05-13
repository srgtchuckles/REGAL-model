---
type: entity
name: REGAL Trial
aliases: [REGAL, GPS-001, NCT03761069]
---

## Overview
REGAL is a Phase 3 randomized controlled trial evaluating galinpepimut-S (GPS) versus best available therapy (BAT) in patients with acute myeloid leukemia (AML) in second complete remission (CR2). It is the pivotal trial for GPS's potential FDA approval and the central subject of this wiki.

## Trial Design

- **Phase:** 3
- **Design:** Open-label, randomized
- **Arms:** GPS vs. BAT (best available therapy)
- **Population:** AML patients in CR2 after successful second-line antileukemic therapy
- **Primary endpoint:** Overall survival (OS)
- **Secondary endpoints:** Not yet populated in wiki
- **Sample size:** ~125–140 patients enrolled; confirmed N=126 (from [[sources/confident-web-7118-regal-dd-feb2026]])
- **Randomization:** ~63 patients per arm (1:1)
- **Sites:** ~95 clinical sites in North America, Europe, and Asia
- **Sponsor:** SELLAS Life Sciences (see [[entities/sellas]])
- **3D Medicines participation:** Up to ~20 additional patients from mainland China possible per Side Letter Agreement (still technically open per 10-Q)

## Statistical Design

- **Interim analysis:** After 60 events (deaths) — conducted
- **Final analysis:** After **80 events (deaths)** — triggers readout
- **IDMC:** Independent Data Monitoring Committee reviews interim + futility
- **Threshold (inferred):** HR < 0.636 at final analysis (from [[sources/confident-web-7118-regal-dd-feb2026]]; not confirmed from protocol documents)
- **Statistical penalty:** None incurred — the one-time event count disclosure in the 10-Q filing does not affect future analyses per SELLAS

## Enrollment Timeline

| Milestone | Date |
|---|---|
| First patient enrolled | ~February 2021 |
| Protocol amendment | ~November 2022 |
| Enrollment complete (126th patient) | April 2024 |

## Event Count Anchors (ALL PUBLIC — FROM SEC FILINGS)

| Date | Event Count | Source |
|---|---|---|
| ~December 2024 | **60** (interim threshold hit) | SELLAS press release |
| January 2025 | IDMC interim analysis completed; "continue without modifications" | SELLAS press release |
| December 26, 2025 | **72** events | [[sources/confident-web-7118-regal-dd-feb2026]] (citing SEC/CRO data) |
| **May 11, 2026** | **78** events | [[sources/sellas-10q-q1-2026]] ← **MOST CURRENT** |

**Final analysis trigger: 80 events.**
**As of May 11, 2026: 2 events remain.**

## Event Velocity Analysis

| Period | Start Events | End Events | Duration | Rate |
|---|---|---|---|---|
| ~Dec 2024 → Dec 26, 2025 | 60 | 72 | ~12 months | 1.0/mo |
| Dec 26, 2025 → May 11, 2026 | 72 | 78 | ~4.5 months | 1.33/mo |

At the recent rate of ~1.33 events/month, the 80th event is expected in approximately **late June to early July 2026** (approximately 6 weeks from the 10-Q filing date). However, with only 2 events needed from ~48 surviving patients, variance is high.

The deceleration from earlier periods (implied >2–3 events/month pre-interim) to the current ~1.33/month is consistent with a cure fraction plateau — the faster-dying non-responders have mostly died; the survivors are dying slowly or not at all.

## IDMC Review History

| Date | Trigger | Outcome |
|---|---|---|
| January 2025 | 60-event interim | "Continue without modifications" |

No additional interim analysis is scheduled. The next analysis is the final (80-event) analysis.

## Current Status (as of May 12, 2026)

- **78/80 events reached** — trial is 2 deaths from final analysis trigger
- Company and investigators remain **fully blinded** to efficacy and survival data
- SELLAS has committed to announcing the 80th event when it occurs
- BLA preparation underway (manufacturing, regulatory consulting)
- No calendar-date trigger provision mentioned in public disclosures

## Regulatory Context

- **GPS designations:** Orphan Drug (AML, MPM, MM — FDA and EMA); Fast Track (AML, MPM, MM — FDA); Rare Pediatric Disease (pediatric AML, Oct 2024 — FDA)
- **BLA pathway:** SELLAS expects REGAL data to form the basis for a BLA submission, subject to statistically significant and clinically meaningful outcome and FDA agreement

## Open Questions

- Exact statistical design (SAP): pre-specified HR threshold, power assumptions
- BAT arm composition: what treatments are allowed; allo-SCT rate (critical confounder)
- Whether REGAL enrolled HLA-A*02:01 enriched patients
- Any calendar-date trigger provision for final analysis if 80th event is delayed
- Full list of secondary endpoints
- Whether 3DMed will enroll mainland China patients (would affect total N)

## Sources
- [[sources/prior-claude-session-apr2026]] — event anchors, cure fraction model context
- [[sources/confident-web-7118-regal-dd-feb2026]] — 72-event anchor (Dec 26, 2025); N=126 confirmation; IDMC "continue without modification" at both reviews
- [[sources/sellas-10q-q1-2026]] — 78-event anchor (May 11, 2026); BLA prep spend; IDMC January 2025 confirmation
- [[sources/jamy-cicic-2025-regal-design]] — trial design paper (not yet summarized in wiki)
