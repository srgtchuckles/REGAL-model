---
type: source
title: "Mapping of novel peptides of WT-1 and presenting HLA alleles that induce epitope-specific HLA-restricted T cells with cytotoxic activity against WT-1+ leukemias (Doubrovina 2012)"
source_type: paper
date_ingested: 2026-04-11
original_file: raw/papers/ (not yet downloaded — citation only)
citation: "Doubrovina E, et al. Blood. 2012 Aug 23;120(8):1633-46. doi: 10.1182/blood-2011-11-394619. PMID: 22623625; PMCID: PMC3429306."
tags: [wt1, hla, cd8, cd4, epitope-mapping, population-coverage]
---

## Summary

Systematic mapping of novel WT-1 peptide epitopes and their HLA restrictions. Identified 41 previously unreported WT-1 epitopes (36 class I, 5 class II). Critical finding for GPS: multiple HLA alleles beyond A*02:01 can present WT-1 peptides and generate cytotoxic T-cell responses, supporting broader population coverage than HLA-A*02:01 alone.

## Key Claims

### Epitopes Identified
- **41 novel WT-1 epitopes total**: 36 MHC class I (CD8+), 5 MHC class II (CD4+)
- Only 1 epitope (126-134 RMFPNAPYL) had been previously described

### HLA Class I Alleles Presenting WT-1 Peptides
- **HLA-A\*02:01**: multiple epitopes including 126-134 RMFPNAPYL
- **HLA-A\*03**: epitopes mapped
- **HLA-A\*24:02**: epitopes mapped
- **HLA-B\*07:02, B\*35:01, B\*40:01, B\*44:02, B\*35:03**: all present WT-1 peptides
- **10 epitopes can be presented by 2–4 different HLA alleles** — broadening coverage

### HLA Class II Alleles Presenting WT-1 Peptides (CD4+)
- **HLA-DRB1\*01:01, DRB1\*04:01, DRB1\*04:02** — CD4+ T-cell responses demonstrated

### Cytotoxic Activity
- T cells generated against 98% of class I epitopes showed HLA-restricted cytotoxicity against peptide-loaded targets
- **75% of 36 evaluated T-cell lines lysed WT-1+ leukemic targets** sharing the restricting HLA allele
- CD4+ T cells from 3 of 5 class II epitopes also showed cytotoxic activity

## Relevance to REGAL Outcome Prediction

### Broader Population Coverage Than Assumed
Our cure fraction model assumed coverage mainly via HLA-A*02:01 (CD8+, ~45% of White AML patients). This paper shows that GPS-style WT-1 targeting could generate cytotoxic responses in patients with HLA-A*03, A*24, multiple B alleles, and HLA-DR alleles.

**Implication:** The HLA-A*02:01-negative component of our cure fraction model (currently 7.3%, based on CD4+ response alone) may be **underestimated**. Patients with HLA-A*03, A*24, or B-allele coverage who mount CD8+ responses against non-A*02 epitopes are not captured in our model.

**However:** GPS as formulated contains specific peptides targeting A*02:01 (the heteroclitic WT1-A1 short peptide). The long peptides (331, 427, 122A1) are intended for CD4+ responses and broader coverage. Whether the GPS formulation covers all the B-allele and A*03/A*24 epitopes identified here requires checking the GPS peptide sequences against this mapping.

### CD4+ Epitopes Mapped to Specific HLA-DR Alleles
HLA-DRB1\*01:01, DRB1\*04:01, DRB1\*04:02 — these are common European alleles. This supports the biological plausibility of the CD4+ response pathway in REGAL's predominantly White patient population.

## Contradictions / Tensions

- The GPS formulation contains 4 specific peptides; whether they include epitopes for A*03, A*24, or B alleles is unclear from this paper alone. If GPS only covers A*02:01 for CD8+, the non-A*02 cytotoxic pathway identified here would not be activated by GPS.
- See [[entities/gps]] for GPS peptide composition details.

## Open Questions

- Do the GPS long peptides (331, 427, 122A1) generate CD8+ responses via HLA-A*03 or A*24 restriction?
- What fraction of the REGAL population would be covered by A*03 + A*24 + B-alleles in addition to A*02:01?

## Sources
- PMC full text fetched: PMC3429306
