---
type: analysis
title: "GPS Cure Fraction Base Case Model"
date: 2026-04-11
sources_used: [prior-claude-session-apr2026]
---

## Question

What is the plausible cure fraction for GPS in REGAL, and what does that imply for the likely hazard ratio and trial outcome?

## Method

Parametric mixture cure model. See [[concepts/mixture-cure-model]] for full methodology.

Cure fraction derived from a chain of Phase 2 data points:

```
π = (HLA-A*02:01+ rate) × (CD8+ response rate) × (DFS plateau in responders)
  + HLA-A*02:01- residual
  + CD4+ durability contribution

= 0.45 × 0.86 × 0.575
+ 0.55 × 0.04
+ 0.018

≈ 26.2%
```

S(t) formula: `S(t) = π + (1 − π) × S_susceptible(t)`

BAT arm modeled as pure Weibull (no cure fraction).

## Findings

| Output | Base Case |
|--------|-----------|
| Aggregate cure fraction (π) | 26.2% |
| Hazard ratio (GPS vs BAT) | ~0.325 |
| GPS 36-month OS | ~42.8% |
| BAT 36-month OS | ~5.6% |
| GPS KM plateau (long-run) | ~26.2% |
| BAT KM plateau | ~0% |

Under base case assumptions, the HR of ~0.325 implies a very large treatment effect.

### Sensitivity: HR as function of π and BAT median OS

The most important sensitivity is BAT median OS (the bear case for GPS):
- BAT mOS = 8 mo (optimistic for GPS): HR ~0.25 → very significant
- BAT mOS = 12 mo (central): HR ~0.325
- BAT mOS = 16 mo (pessimistic for GPS): HR moves toward 0.5+ depending on π

At π = 10% (conservative cure fraction), even with BAT mOS = 16mo, the model still suggests HR < 0.75 — but the margin narrows considerably.

## Confidence Level

**Low–Medium** — the directional signal (meaningful cure fraction exists) is supported by Phase 2 biology and curve shape. The specific numbers (26.2%, HR 0.325) are **not reliable point estimates** because:
- Phase 2 n is 9 responders vs 5 non-responders (underpowered, p>0.05)
- Figure 5 DFS plateau has wide uncertainty at 65+ months with few patients remaining
- Cure fraction estimate has a known double-counting bias (upward)
- BAT arm performance is unknown and could be materially higher than assumed

## Implications for Outcome Prediction

If Phase 2 signal is real and scales to REGAL's larger n:
- GPS almost certainly meets its primary OS endpoint
- HR in the range 0.3–0.5 would be a landmark result in AML CR2

If Phase 2 signal is noise (small n artifact):
- GPS cure fraction is near zero, trial looks like two Weibull curves with modest or no separation
- Trial fails to meet primary endpoint

The central question REGAL will answer: does the Phase 2 immune responder survival advantage reflect a true biological signal or small-n selection artifact?

## Caveats

1. No actual REGAL KM data available — model is built entirely from Phase 2 (n=22) and public event count anchors
2. BAT arm allo-SCT rate is unknown — a high allo-SCT rate could substantially raise BAT arm OS
3. HLA enrichment of REGAL enrollment unknown — if REGAL enrolled HLA-A\*02:01+ enriched patients, π estimate should be revised upward
4. Python model exists at `~/Claude_REGAL_model.py` — needs path fix for local output confirmed working as of Apr 5 2026

## Related Pages

- [[concepts/mixture-cure-model]]
- [[entities/gps]]
- [[entities/regal-trial]]
- [[entities/bat-arm]]
- [[entities/aml-cr2]]
