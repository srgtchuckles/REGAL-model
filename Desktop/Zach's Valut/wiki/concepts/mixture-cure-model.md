---
type: concept
name: Parametric Mixture Cure Model
aliases: [cure fraction model, split population model, MCM]
---

## Definition

A parametric mixture cure model (MCM) assumes the trial population is composed of two latent subgroups:
- A fraction **π** ("cured") — hazard is zero, they will never experience the event
- A fraction **(1 − π)** ("susceptible") — will eventually have the event; survival follows a parametric distribution

The overall survival function is derived from the law of total probability:

```
S(t) = π + (1 − π) × S_susceptible(t)
```

As t → ∞, S_susceptible(t) → 0, so the curve plateaus at π. This plateau is the visual signature of a cure fraction in a KM curve.

## How It Applies to REGAL

GPS's mechanism (WT1-specific T-cell immunosurveillance) provides biological justification for a cure fraction: patients who mount durable immune responses may achieve persistent remission maintenance with effectively zero ongoing hazard. The Phase 2 Figure 5 DFS curve for immune responders shows exactly this plateau signature (~55–60% at 65+ months).

The model structure for REGAL:
- **GPS arm:** Mixture cure model with π ~26.2% (base case), Weibull susceptible component (median OS ~20mo, shape k~1.3)
- **BAT arm:** Pure Weibull (no cure fraction modeled; possible underestimate if allo-SCT is common in BAT)

## Weibull Choice for S_susceptible

Weibull is used because:
- Two parameters (scale λ, shape k) can be fit from median OS + one shape assumption
- k > 1 gives increasing hazard over time (appropriate for AML relapse biology)
- Reduces to exponential when k = 1

## Key Parameters and Their Sources

| Parameter | Base Case | Source |
|-----------|-----------|--------|
| π (cure fraction) | 26.2% | Derived from Phase 2 HLA + immune response chain |
| HLA-A\*02:01+ rate | 45% | Hardy-Weinberg from allele frequency literature |
| CD8+ response rate | 86% | Phase 2 n=7 HLA+ patients |
| DFS plateau (responders) | 57.5% | Figure 5, Phase 2 paper |
| GPS non-cured median OS | ~20 mo | Assumed; needs grounding from data |
| BAT median OS | 8–16 mo | Historical range; central estimate ~12 mo |

## Sensitivity Analysis Dimensions

The two most important parameters to stress-test:
1. **π** (cure fraction): range 10–25% covers most plausible scenarios
2. **BAT median OS**: range 8–16 months; the higher end (BAT performing well) is the bear case for GPS

A 2D HR sensitivity heatmap sweeping these two parameters (implemented in the Python model) is the key analytical output.

## Limitations and Caveats

- Phase 2 sample: n=9 responders vs n=5 non-responders — extremely small; p-values non-significant (DFS p=0.11, OS p=0.08)
- Double-counting issue: Figure 5 plateau is from a mixed HLA immune responder population; applying it to only the CD8+ chain slightly inflates the aggregate cure fraction
- BAT arm has no cure fraction modeled — may understate BAT performance if allo-SCT is common
- Two-anchor fitting (60 events at mo 46, 72 at mo 58) is insufficient for rigorous model validation; digitized KM curves are needed

## Related Pages

- [[sources/prior-claude-session-apr2026]]
- [[entities/gps]]
- [[entities/regal-trial]]
- [[entities/bat-arm]]
- [[analyses/cure-fraction-base-case]]
