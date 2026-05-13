---
type: source
title: "WT1 peptide vaccination in AML induces short-lived WT1-specific immune responses (Uttenthal 2014)"
source_type: paper
date_ingested: 2026-04-11
original_file: raw/papers/ (not yet downloaded — citation only)
citation: "Uttenthal B, et al. Br J Haematol. 2014 Feb;164(3):366-75. doi: 10.1111/bjh.12637. PMID: 24422723; PMCID: PMC4253125."
tags: [wt1, vaccine, immune-response, short-lived, cd8, aml, negative-signal]
---

## Summary

A WT1 peptide vaccine study showing **short-lived** immune responses and **no clinical benefit** in AML patients. This is an important negative/cautionary signal. Critically, the vaccine used was **different from GPS** (only 2 peptides, different adjuvant), and the authors propose specific mechanistic reasons why responses were transient — reasons that GPS's design partially addresses.

## Key Claims

### Patient Population
- 8 HLA-A2-positive patients with poor-risk AML (CR1, CR2, PR, or progressive disease)
- Median age 65 years

### Vaccine Used (NOT GPS)
- **pWT126 and pWT235** (2 HLA-A2-binding WT1 epitopes) + PADRE (pan-DR helper epitope)
- Adjuvant: **Montanide** (no GM-CSF)
- 5 vaccination cycles at 3-weekly intervals, dose escalation (0.3, 0.6, 1.0 mg)
- **Key difference from GPS**: only 2 peptides, no GM-CSF, different adjuvant

### Immune Response Rates
- **86% response rate**: 6/7 evaluable patients showed WT1-specific CD8+ T cells (tetramer + ELISPOT)
- Responses detected within 2 weeks of initial vaccination
- **Critical: responses were SHORT-LIVED and did not expand upon re-stimulation**

### Clinical Outcomes
- **No correlation** between immune response, WT1 mRNA reduction, or clinical benefit
- Only 2/5 evaluable patients showed WT1 mRNA reduction coinciding with immune responses

### Mechanistic Explanation for Short-Lived Responses
Authors propose 3 factors:
1. **Adjuvant effect**: Sustained peptide presentation without continuous helper/danger signal may induce CD8+ T-cell tolerance rather than memory
2. **Low-avidity T cells**: High-avidity clones may be deleted due to WT1's expression in normal tissues during development (central tolerance)
3. **Limited avidity maturation**: Self-antigens like WT1 restrict memory CD8+ T-cell avidity maturation vs. foreign antigens

## Relevance to REGAL Outcome Prediction

### Important Negative Signal — But Applies to a Different Vaccine
This paper shows that a 2-peptide WT1 vaccine with Montanide-only adjuvant generates transient, functionally impaired CD8+ responses with no clinical benefit. This is the key cautionary counterpoint to our optimistic GPS model.

**However, GPS is specifically designed to overcome these limitations:**
- **GM-CSF adjuvant**: GPS uses GM-CSF in addition to Montanide — GM-CSF provides the "danger signal" that the Uttenthal paper identifies as missing. This may be why GPS Phase 2 showed durable responses while this vaccine did not.
- **4 peptides vs 2**: The multi-peptide design of GPS (including 2 native long peptides for CD4+ help and 1 heteroclitic long peptide for both CD4+/CD8+) creates the sustained helper T-cell signal that addresses limitation #1.
- **CD4+ help**: GPS explicitly targets CD4+ helper T cells via long peptides. The Maslak Phase 2 data shows that CD4+ responders never relapsed — GPS's CD4+ arm may be the key mechanism preventing the tolerance/exhaustion seen in this paper.

### Avidity Maturation Concern
The concern about low-avidity T-cell selection (WT1 as self-antigen) is a real biological risk. The heteroclitic WT1-A1 peptide in GPS was specifically designed to generate higher-avidity responses than the native peptide — this addresses limitation #2 directly. See [[sources/zirlik-2006-heteroclitic]].

### Downside Risk for Model
If GPS's GM-CSF + multi-peptide design is insufficient to overcome these tolerance mechanisms, the Phase 2 DFS plateau may erode in REGAL's larger, more heterogeneous population. The Phase 2 patients may have been a selected, immunologically favorable group.

## Contradictions / Tensions

Directly contradicts the Phase 2 GPS durability signal. Reconciliation requires believing that GPS's specific design (GM-CSF + heteroclitic peptide + CD4+ epitopes) genuinely overcomes the tolerance mechanisms identified here. This is plausible but not proven.

## Open Questions

- Are there head-to-head comparisons of GPS vs. similar 2-peptide + Montanide vaccines in matched patient populations?
- Do GPS Phase 2 responders show avidity maturation (high-avidity T cells) vs. the low-avidity cells in this study?

## Sources
- PMC full text fetched: PMC4253125
