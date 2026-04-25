---
type: entity
name: Galinpepimut-S (GPS)
aliases: [GPS, galinpepimut-S, WT1 vaccine, SLS009]
---

## Overview
Galinpepimut-S (GPS) is a WT1 antigen-targeting immunotherapy (multi-peptide cancer vaccine) developed by Sellas Life Sciences. It targets Wilms' Tumor protein 1 (WT1), which is overexpressed in AML blasts. GPS is designed to stimulate both CD4+ and CD8+ T-cell responses against WT1-expressing tumor cells.

## Key Facts

- Drug class: Multi-peptide cancer vaccine / immunotherapy
- Target: WT1 (Wilms' Tumor 1 protein) — overexpressed on AML blasts and leukemic stem cells
- Administration: Subcutaneous injection
- Prior development: Memorial Sloan Kettering; acquired by Sellas Life Sciences

## Mechanism of Action

GPS presents WT1-derived peptides to prime both CD8+ cytotoxic T-cells (via MHC Class I) and CD4+ helper T-cells (via MHC Class II). The CD8+ arm directly kills WT1-expressing leukemic blasts; the CD4+ arm sustains long-term immune memory and surveillance. This dual activation is why a durable cure fraction is biologically plausible.

**Key mechanistic detail:** GPS uses a *heteroclitic* peptide (WT1-A1) that generates immune responses cross-reactive to the native WT1-A peptide present in leukemic cells. Cross-reactivity confirmed in all 5 ELISPOT-tested patients in Phase 2.

**HLA restriction:** CD8+ responses require HLA-A*02:01 (or compatible allele) for MHC Class I peptide presentation. ~45% of predominantly White AML populations carry HLA-A*02:01. CD4+ responses require HLA-DR alleles (less well characterized for GPS).

## Phase 2 Immunologic Correlate Data

*From Phase 2 publication (n=22 enrolled, n=14 with immunologic data):*

- **CD8+ response rate:** 6/7 HLA-A\*02:01+ patients tested = **86%**
- **CD4+ response rate:** 4/9 tested = **44%**
  - None of the CD4+ responders relapsed (most striking finding in the dataset)
- **DFS by immune response (Figure 5, Panel A):**
  - Responders (n=9): median not reached, plateau ~55–60% at 65+ months
  - Non-responders (n=5): median 15.6 months, tail ~20%
  - P=0.11 (non-significant — underpowered)
- **OS by immune response (Figure 5, Panel B):**
  - Responders: median not reached, plateau ~65–70% at 65+ months
  - Non-responders: median 35.8 months, tail ~20–25%
  - P=0.08 (non-significant — underpowered)

## Clinical History
*To be expanded from additional sources.*

## Safety Profile
*To be populated from sources.*

## Relevance to REGAL Outcome Prediction
The Phase 2 immunologic data is the foundational input for the cure fraction model. The key question REGAL answers: is the Phase 2 immune responder survival advantage (DFS plateau ~57.5%) a real biological signal or a small-n artifact? See [[analyses/cure-fraction-base-case]].

## Open Questions
- Was REGAL enrollment enriched for HLA-A\*02:01+ patients? (Critical for cure fraction estimate)
- Are there published KM curves from Phase 1 GPS data that can be digitized?
- What specific HLA-DR alleles are the GPS CD4+ epitopes designed around?
- Does CD4+ response remain the strongest predictor of durable remission in larger series?

## Sources
- [[sources/prior-claude-session-apr2026]]
