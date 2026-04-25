---
type: source
title: "Prior Claude Session — REGAL/GPS Cure Fraction Analysis (Apr 4–5 2026)"
source_type: chat
date_ingested: 2026-04-11
original_file: raw/other/claude-chat-REGAL-trial-code.rtfd
tags: [gps, regal, aml-cr2, cure-fraction, mixture-cure-model, hla, phase2, statistics]
---

## Summary

This was an extended analysis session building a parametric mixture cure model for the REGAL trial outcome, grounded in Phase 2 GPS immunologic correlate data and HLA population frequency literature. The session progressed from critiquing a Reddit cure-fraction argument, to extracting Phase 2 immunologic data, to building a full Python model with sensitivity analysis.

A key early error: Claude initially misidentified REGAL as the imetelstat/MDS trial; the user corrected this. The rest of the session is specific to GPS/AML CR2.

## Key Claims

### Public Event Anchors (from SEC filings)
- 60 events at month 46
- 72 events at month 58
- Event velocity deceleration between these anchors is observable and was the basis of the Reddit post

### Phase 2 Immunologic Data (n=22 enrolled; n=14 with immunologic data)
- **CD8+ response (HLA-A\*02:01+ patients):** 6/7 tested = **86% response rate**
  - Cross-reactivity to native WT1-A confirmed in all 5 ELISPOT-tested patients (critical: heteroclitic peptide generates response to native leukemic antigen)
  - 4/5 tetramer-tested patients showed increased WT1-A tetramer-positive CD8+ cells post-vaccination
- **CD4+ response:** 4/9 tested = **44% response rate**
  - No HLA-DR subtype correlation detectable (n too small)
  - **Most important finding: none of the CD4+ responders relapsed**
- **DFS/OS by immune response (Figure 5):**
  - Responders (n=9): DFS median not reached, plateau ~55–60% at 65+ months; OS median not reached, plateau ~65–70%
  - Non-responders (n=5): DFS median 15.6 months, tail ~20%; OS median 35.8 months, tail ~20–25%
  - P-values: DFS p=0.11, OS p=0.08 — curve separation visually dramatic but statistically non-significant (underpowered)

### HLA-A\*02:01 Population Frequencies
- European/White: ~27% allele frequency → ~45% carrier rate (Hardy-Weinberg)
- African/Black: ~12% allele frequency → ~22% carrier rate
- Asian/Pacific Islander: ~6.5% allele frequency → ~13% carrier rate
- REGAL enrollment skews White/elderly (AML demographics), so ~45% carrier rate is the relevant estimate

### Derived Cure Fraction (Base Case: 26.2%)
Chain of reasoning:
1. ~45% of trial is HLA-A\*02:01+ → eligible for strong CD8+ response
2. × 86% mount a CD8+ immune response
3. × 57.5% (midpoint of 55–60% DFS plateau among responders) = **22.2% from CD8+ chain**
4. HLA-A\*02:01− residual contribution: 0.55 × ~4% = **~2.2%**
5. CD4+ durability bonus: **~1.8%**
6. Total: **~26.2%**

**Known limitation:** The DFS plateau in Figure 5 is from a mixed HLA population of immune responders, so applying it only to the CD8+ chain slightly double-counts. Biases cure fraction upward by a few percentage points.

### Model Formula
Parametric mixture cure model (split population model):

```
S(t) = π + (1 − π) × S_susceptible(t)
```

- π = cure fraction (~26.2% base case)
- S_susceptible(t) = Weibull survival for non-cured patients
- As t → ∞, S_susceptible → 0, so curve plateaus at π
- Weibull shape k > 1 (increasing hazard appropriate for AML)

### Base Case Model Outputs
- Aggregate cure fraction: 26.2%
- Hazard ratio: **~0.325** (strongly favoring GPS under base case)
- GPS 36-month OS: ~42.8%
- BAT 36-month OS: ~5.6%
- BAT has no cure fraction modeled (pure Weibull, hits ~0% by month 48)

### Python Model
A working Python implementation exists at `~/Claude_REGAL_model.py` on the user's machine. Features: 6-panel figure (survival curves, hazard, event velocity, derivation box, HR heatmap, p-value sensitivity grid). Fixed save path issue (was using `/mnt/user-data/outputs/` container path).

### Reddit Post Critique (r/pennystocks, user Confident-Web-7118)
- Core observation valid: event deceleration pattern is real and worth noting
- "99.99% probability of success" claim: not defensible — rhetorical, not statistical
- Fitting 6 model families to 2 anchors is overdetermined (curve fitting, not model validation)
- Biological cap argument (AML CR2 has 6–12 month cap, so long-tail survivors must be cured): overstated — frailty, patient selection, and supportive care variability can produce apparent tails without true cure fraction
- Ignored control arm entirely — cure fraction in treatment arm only means nothing without comparative KM shape

## Relevance to REGAL Outcome Prediction

This session built the foundational quantitative model for the prediction. The cure fraction estimate (~26.2% base case, range ~10–25% stress-tested) is the central parameter. The HR of ~0.325 under base case assumptions suggests a strongly positive trial if the Phase 2 signal scales — but the key uncertainty is whether the tiny Phase 2 n (9 vs 5) represents a real signal or noise.

## Contradictions / Tensions

- The model's base case HR (0.325) implies an extremely strong treatment effect — unusually large for a Phase 3 oncology trial. This should prompt skepticism: either GPS is genuinely exceptional, or the Phase 2 data is upward-biased due to small n and potential patient selection.
- No BAT arm cure fraction was modeled, but in practice some BAT patients may proceed to allo-SCT (potentially curative). This could materially improve the BAT arm. See [[entities/bat-arm]].

## Open Questions

- Was REGAL enrollment enriched for HLA-A\*02:01+ patients? (Would shift aggregate cure fraction estimate up)
- Do any conference abstracts hint at event velocity or interim DSMB activity?
- Are there digitized KM curves from published GPS publications for model fitting beyond two event anchors?
- What is the actual sample size and power assumption in the REGAL statistical analysis plan?
- What specific therapies were used in the BAT arm and at what frequency?
