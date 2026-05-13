---
type: source
title: "REGAL: galinpepimut-S vs. best available therapy as maintenance therapy for AML in second remission (Jamy & Cicic 2025)"
source_type: paper
date_ingested: 2026-04-11
original_file: raw/papers/ (not yet downloaded — citation only)
citation: "Jamy O, Cicic D. Future Oncol. 2025 Jan;21(1):73-81. doi: 10.1080/14796694.2024.2433935. PMID: 39606837; PMCID: PMC11760237."
tags: [regal, gps, trial-design, aml-cr2, bat-arm, endpoints, statistics]
---

## Summary

This is a Future Oncology paper describing the REGAL Phase 3 trial design in detail. It is the single most important source for understanding the trial mechanics. Key findings that materially affect the outcome prediction model:

1. **Allo-SCT is excluded from the BAT arm** — patients with planned allo-SCT are ineligible for REGAL. This eliminates the biggest model uncertainty (that BAT arm patients could receive a potentially curative transplant).
2. **Protocol's own BAT median OS assumption is 8.0 months** — our model used 12.0mo, which was actually conservative (pessimistic for GPS). The trial's own power calculation assumes a more favorable GPS advantage.
3. **Interim analysis at 60 deaths** — matches the SEC-filed event anchor exactly. The interim analysis has already occurred.
4. **Target final analysis: 80 deaths** — only 8 events beyond the 72 reported at mo58 are needed for the final readout.

## Key Claims

### Primary Endpoint
- Overall survival (OS) in AML CR2/CRp2
- Cox PH model, stratified by randomization factors
- 1-sided α = 0.025; H0: HR ≥ 1 vs H1: HR < 1

### Secondary Endpoints
- Leukemia-Free Survival (LFS)
- OS and LFS rates at 6, 9, 12 months
- MRD assessment (multiple methods)
- Safety and tolerability

### Sample Size & Power
- **~125–140 patients enrolled, 1:1 randomization**
- **80 deaths required for final analysis**
- **90% power** assuming HR = 0.636, BAT mOS = 8.0 months, GPS mOS = 12.6 months
- Protocol's own HR assumption for powering: **0.636**

### Interim Analysis
- **1 pre-specified efficacy interim analysis after 60 deaths** (Lan-DeMets O'Brien-Fleming alpha spending)
- IDMC makes stopping recommendations
- **The interim analysis at 60 deaths matches the SEC-filed anchor (60 events at mo46)** — the interim has occurred

### BAT Arm Definition (CRITICAL)
Permitted BAT therapies:
- Observation (including hydroxyurea)
- Hypomethylating agent: decitabine or azacitidine, alone or with venetoclax
- Low-dose cytarabine (LDAC)

**Excluded / ineligible:**
- Patients with **imminently planned allo-SCT are ineligible for enrollment**
- Patients on molecularly targeted agents (continuation)
- Active CNS leukemia

### GPS Dosing Schedule
- Induction: every 2 weeks × 6 doses (weeks 0–10)
- Booster 1: every 4 weeks × 6 doses (weeks 14–34)
- Booster 2: every 6 weeks × 3 doses (weeks 40–52)

### Stratification Factors
1. CR1 duration (≥12 months vs. <12 months)
2. Baseline cytogenetics risk (poor vs. all other)
3. CR2 vs. CRp2 status
4. MRD presence/absence after second remission

### Key Inclusion Criteria
- Age ≥18, AML per WHO criteria (de novo or secondary)
- CR2/CRp2 for relapsed AML
- ≥300 lymphocytes/µL
- **Not a candidate for allo-SCT**
- Enrolled within 6 months of achieving CR2/CRp2
- ECOG PS 0–3

## Relevance to REGAL Outcome Prediction

### Allo-SCT Exclusion
This is the single most important finding from this paper. Our model's largest uncertainty was whether BAT arm patients could receive allo-SCT (potentially curative). **They cannot — patients eligible for allo-SCT are excluded from enrollment.** The BAT arm is definitionally limited to non-curative therapies. This removes the main bear-case scenario for GPS.

### BAT mOS Assumption
The protocol powers the study assuming BAT mOS = 8.0 months. Our base case model used 12.0 months — we were being conservative. Under the protocol's own assumption (8.0mo), the expected GPS benefit is even larger. The sensitivity heatmap's most pessimistic BAT mOS of 16mo is likely unrealistically high given the protocol's own 8.0mo estimate and the exclusion of allo-SCT.

### Interim Analysis Already Occurred
The interim analysis at 60 events has already happened (60 events at month 46 per SEC filing). The trial was not stopped for futility or overwhelming efficacy. This constrains the cure fraction range: the IDMC continued the trial, which means the interim result was likely in the "continue" zone (neither stop for futility nor stop for overwhelming efficacy). The O'Brien-Fleming boundary at interim typically requires a very extreme result to stop early — continuation is the most common outcome.

### 8 Events Needed for Final Analysis
Only 72 events reported at mo58; 80 needed for final analysis. With event velocity slowing, readout timing depends heavily on how far along the cure fraction plateau has developed.

## Contradictions / Tensions

- The protocol assumes BAT mOS = 8.0mo and GPS mOS = 12.6mo — a much more modest GPS benefit than our model's base case (HR ~0.343). If the protocol's own assumptions are close to reality, the trial succeeds but with a less dramatic effect size than the cure fraction model predicts.
- The HR threshold for 90% power is 0.636, not the 0.325–0.343 from the cure fraction model. This suggests the protocol designers were more conservative about GPS's benefit than our model.
- BAT arm is now restricted to HMAs, VEN, LDAC, or observation — the actual 2020–2025 AML landscape may have made some of these more effective than the 2019 protocol assumption, which could push realized BAT mOS above 8.0mo.

## Open Questions

- Did the IDMC at the 60-event interim cross any pre-specified efficacy boundary, or simply continue?
- What was the actual enrollment period and patient count achieved (125 or 140)?
- Are there conference abstracts (ASH, EHA) with any REGAL interim data or event velocity updates?
- What is the event velocity in the post-mo58 period — are events still accumulating slowly (cure fraction plateau) or accelerating?

## Sources
- PMC full text fetched: PMC11760237
