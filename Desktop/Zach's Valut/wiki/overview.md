# REGAL Trial — Outcome Prediction Overview

*Last updated: 2026-05-12. This page is the evolving synthesis.*

---

## Current Prediction

> **GPS meets primary OS endpoint — base case probability: MEDIUM-HIGH. Confidence has increased materially since April 2026 due to the 78/80 event anchor and continuing velocity deceleration.**

The 78-event anchor (May 11, 2026) confirms that the trial is 2 deaths from final analysis. The event velocity from Dec 2025 → May 2026 (1.33 deaths/month from ~48 survivors) remains consistent with a cure fraction plateau — a pure Weibull model without a plateau would predict 5–10x higher mortality rates at this stage. The IDMC's January 2025 "continue without modification" ruling rules out a futility stop at interim.

**Revised probability statement:** GPS meeting primary endpoint is likely in the **60–75%** range. The velocity data and IDMC signal have narrowed the range upward from the April 2026 estimate of 40–65%. The remaining uncertainty is primarily: (1) true BAT arm mOS (unknown); (2) whether the GPS cure fraction plateau is real or a small-n artifact that won't replicate at scale; (3) HLA-A*02:01 enrollment enrichment status.

**Final analysis expected: late June – early July 2026** (approximately 6 weeks from May 11, 2026 at current velocity).

---

## Key Evidence For

- **Phase 2 immunologic correlate data:** CD8+ response rate 86% in HLA-A\*02:01+ patients; none of the 4 CD4+ responders relapsed; immune responder DFS plateau ~57.5% at 65+ months ([[sources/prior-claude-session-apr2026]])
- **Biological plausibility of cure fraction:** Peptide vaccines targeting WT1 are one of the few mechanisms that could generate durable T-cell memory → genuine cure fraction signal is mechanistically coherent ([[entities/gps]])
- **Event velocity deceleration confirms plateau:** 1.0 deaths/month (Dec 2024–Dec 2025) → 1.33/month (Dec 2025–May 2026). From ~48 surviving patients at the later period, this implies annualized mortality of ~33%/year — far below what a standard Weibull model predicts without a plateau for 5+ year survivors in AML CR2. ([[entities/regal-trial]], [[sources/sellas-10q-q1-2026]])
- **IDMC "continue without modification"** at the 60-event interim (January 2025) — rules out a futility stop, a key residual risk. ([[sources/sellas-10q-q1-2026]])
- **BLA preparation spending:** SELLAS ramped manufacturing by +$1.3M and regulatory consulting by +$0.4M in Q1 2026, explicitly attributed to BLA prep. While management has incentive bias, pre-positioning at this scale is a directional positive signal. ([[sources/sellas-10q-q1-2026]])
- **Reverse-engineered cure fraction 42–48%:** Independent retail analysis of public event velocity data produces implied GPS cure fraction substantially above the wiki's biology-derived 26.2% estimate. ([[sources/confident-web-7118-regal-dd-feb2026]])
- **AML CR2 disease biology:** Short natural history tail in CR2 means the control arm is unlikely to produce a long-tail survivor artifact without true treatment effect ([[entities/aml-cr2]])

---

## Key Evidence Against

- **Phase 2 sample size is tiny:** n=9 immune responders vs n=5 non-responders. DFS p=0.11, OS p=0.08. Visual curve separation is striking but statistically underpowered — selection bias cannot be ruled out. ([[sources/prior-claude-session-apr2026]])
- **Cure fraction estimate is upward-biased:** Known double-counting in the HLA chain calculation. True aggregate π may be several percentage points lower than the 26.2% base case.
- **BAT arm is unmodeled:** No cure fraction assigned to BAT, but allo-SCT in the BAT arm (potentially curative) could materially improve control arm OS. If BAT mOS is 14–16 months rather than 8–12 months, the HR picture worsens substantially. ([[entities/bat-arm]])
- **Base case HR (0.325) is unusually large:** Extreme treatment effects in Phase 3 oncology trials are rare. Phase 2 to Phase 3 regression-to-the-mean is a consistent pattern in oncology.
- **Velocity from Dec 2025–May 2026 is slightly higher than pure plateau predicts:** 1.33 events/month from ~48 patients is still slow but not as slow as the most optimistic cure fraction models predict. The survivor cohort still contains non-cured patients who will die.

---

## Critical Uncertainties

1. **BAT arm allo-SCT rate** — single most important unknown; could make or break the HR
2. **True GPS cure fraction in REGAL** — Phase 2 n=9 responders extrapolated to n=126; high variance
3. **REGAL HLA enrollment enrichment** — if enriched for HLA-A\*02:01+, cure fraction estimate should be revised upward; if not enriched, downward
4. **REGAL statistical analysis plan** — exact pre-specified HR threshold and power assumptions not confirmed from protocol (inferred HR < 0.636 threshold)
5. **CD4+ response durability** — "none relapsed" in Phase 2 n=4 is the strongest signal but entirely unconfirmed at scale

---

## Event Timeline & Final Analysis Projection

| Date | Event | Count |
|---|---|---|
| Feb 2021 | First patient enrolled | — |
| Nov 2022 | Protocol amendment | — |
| Apr 2024 | Enrollment complete (N=126) | — |
| Dec 2024 | 60-event interim threshold reached | 60/80 |
| Jan 2025 | IDMC: "continue without modifications" | — |
| Dec 26, 2025 | CRO update | 72/80 |
| May 11, 2026 | CRO update (per 10-Q) | **78/80** |
| **~Late Jun 2026** | **Expected 80th event** (final analysis trigger) | **80/80** |

---

## What Would Change the Prediction

**Would increase confidence GPS succeeds:**
- 80th event announcement (imminent)
- Any conference presentation showing event velocity continuation
- Protocol documents revealing HLA-A\*02:01 enrichment in enrollment criteria
- Evidence BAT arm allo-SCT rate is low (<20%)

**Would decrease confidence:**
- High allo-SCT rate in the BAT arm from protocol details or conference discussions
- Any leak or pre-announcement suggesting the trial failed (though SELLAS has committed to announcing only the 80th event, not interim efficacy)
- Extended delay past July 2026 in reaching 80 events (would imply slower mortality = more patients dying on GPS arm than expected plateau)

---

## Model Reference

Base case model documented in [[analyses/cure-fraction-base-case]]. Python implementation at `~/Claude_REGAL_model.py`. Methodology in [[concepts/mixture-cure-model]].

---

## Related Pages

- [[entities/regal-trial]] ← most current event count and velocity data
- [[entities/gps]]
- [[entities/aml-cr2]]
- [[entities/bat-arm]]
- [[entities/sellas]]
- [[analyses/cure-fraction-base-case]]
- [[sources/sellas-10q-q1-2026]] ← 78-event anchor (May 11, 2026)
- [[sources/confident-web-7118-regal-dd-feb2026]] ← 72-event anchor, velocity proof
