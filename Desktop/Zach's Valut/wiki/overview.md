# REGAL Trial — Outcome Prediction Overview

*Last updated: 2026-04-11. This page is the evolving synthesis.*

---

## Current Prediction

> **GPS meets primary OS endpoint — base case probability: MEDIUM-HIGH under biologically optimistic assumptions, but confidence is LOW due to small Phase 2 n.**

The mixture cure model built from Phase 2 data produces a base case HR of ~0.325 (GPS vs BAT), implying a very large treatment effect. If this signal scales to REGAL's larger n, the trial succeeds convincingly. The critical uncertainty is whether the Phase 2 immune responder survival plateau (~57.5% DFS at 65+ months, n=9) is a real biological signal or a small-n artifact.

**Probability statement:** GPS meeting primary endpoint is plausible — perhaps 40–65% — but this range is wide because the Phase 2 evidence is directionally strong and biologically coherent, yet statistically underpowered. A proper confidence interval requires more sources.

---

## Key Evidence For

- **Phase 2 immunologic correlate data:** CD8+ response rate 86% in HLA-A\*02:01+ patients; none of the 4 CD4+ responders relapsed; immune responder DFS plateau ~57.5% at 65+ months ([[sources/prior-claude-session-apr2026]])
- **Biological plausibility of cure fraction:** Peptide vaccines targeting tumor-specific antigens (WT1) are one of the few mechanisms that could generate durable T-cell memory → genuine cure fraction signal is mechanistically coherent ([[entities/gps]])
- **AML CR2 disease biology:** Short natural history tail in CR2 means the control arm is unlikely to produce a long-tail survivor artifact without true treatment effect — strengthening the interpretation of the GPS responder plateau ([[entities/aml-cr2]])
- **Event deceleration observed:** The 60→72 event trajectory from SEC filings shows velocity deceleration consistent with a cure fraction plateau developing ([[entities/regal-trial]])

---

## Key Evidence Against

- **Phase 2 sample size is tiny:** n=9 immune responders vs n=5 non-responders. DFS p=0.11, OS p=0.08. The visual curve separation is striking but statistically underpowered — selection bias cannot be ruled out.
- **Cure fraction estimate is upward-biased:** Known double-counting in the HLA chain calculation. True aggregate cure fraction may be several percentage points lower than the 26.2% base case.
- **BAT arm is unmodeled:** No cure fraction assigned to BAT, but allo-SCT in the BAT arm (potentially curative) could materially improve control arm OS. If BAT mOS is 14–16 months rather than 8–12 months, the HR picture worsens substantially.
- **Base case HR (0.325) is unusually large:** Extreme treatment effects in Phase 3 oncology trials are rare. This should prompt skepticism about whether Phase 2 numbers will hold.

---

## Critical Uncertainties

1. **BAT arm allo-SCT rate** — single most important unknown; could make or break the HR
2. **REGAL HLA enrollment enrichment** — if enriched for HLA-A\*02:01+, cure fraction estimate should be revised upward
3. **True Phase 2 DFS plateau** — only n=9 responders with long follow-up; the 57.5% figure has wide CIs
4. **REGAL statistical analysis plan** — sample size, power assumptions, pre-specified HR threshold, interim analysis rules
5. **CD4+ response durability** — the "none relapsed" finding is the strongest signal in Phase 2 but entirely unconfirmed with only n=4

---

## What Would Change the Prediction

**Would increase confidence GPS succeeds:**
- Evidence that REGAL enrolled HLA-A\*02:01+ enriched patients
- DSMB continuation signals (no interim futility stop)
- Conference abstracts showing event deceleration continuing beyond month 58
- Published Phase 1 KM curves with similar immune responder plateau

**Would decrease confidence:**
- High allo-SCT rate in the BAT arm (from protocol details or conference discussions)
- Any futility interim analysis that was triggered
- New papers showing GPS CD8+ response rates lower than 86% in broader populations
- Evidence that the Phase 2 immune responder group was selected on disease biology rather than vaccine effect

---

## Model Reference

Base case model documented in [[analyses/cure-fraction-base-case]]. Python implementation at `~/Claude_REGAL_model.py`. Methodology in [[concepts/mixture-cure-model]].

---

## Related Pages

- [[entities/regal-trial]]
- [[entities/gps]]
- [[entities/aml-cr2]]
- [[entities/bat-arm]]
- [[analyses/cure-fraction-base-case]]
