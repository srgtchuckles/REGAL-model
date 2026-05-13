# Wiki Log

Append-only activity log. Format: `## [YYYY-MM-DD] operation | description`

---

## [2026-04-11] init | Wiki scaffolded for REGAL trial analysis
Created directory structure, CLAUDE.md schema, index.md, log.md, overview.md, and 5 entity stubs (GPS, REGAL trial, AML CR2, BAT arm, Sellas).

## [2026-04-11] ingest | Prior Claude session (Apr 4–5 2026) — cure fraction model
Ingested RTF export of prior Claude.ai session. Extracted: Phase 2 immunologic correlate data (CD8+ 86%, CD4+ 44%, Figure 5 DFS plateau ~57.5%), HLA-A*02:01 population frequencies, mixture cure model derivation (π~26.2%, HR~0.325), event anchors from SEC filings (60@mo46, 72@mo58). Created: sources/prior-claude-session-apr2026.md, concepts/mixture-cure-model.md, analyses/cure-fraction-base-case.md. Updated: entities/gps.md, entities/regal-trial.md, overview.md, index.md.

## [2026-05-12] ingest | SELLAS 10-Q Q1 2026 — 78/80 events; BLA prep; financials
Ingested SELLAS Life Sciences Form 10-Q for quarter ended March 31, 2026 (filed May 12, 2026). CRITICAL: CRO confirmed 78 events as of May 11, 2026 (final analysis triggers at 80). Final readout expected ~late June–July 2026. Event velocity Dec 2025–May 2026: 1.33/month from ~48 survivors, consistent with cure fraction plateau. BLA prep spend: +$1.3M manufacturing, +$0.4M regulatory consulting in Q1 2026 (explicitly attributed to REGAL BLA prep). Cash: ~$114.6M. Created: sources/sellas-10q-q1-2026.md (+ PDF copied to wiki/sources/). Updated: entities/sellas.md (fully populated), entities/regal-trial.md (event table, velocity analysis), overview.md (probability revised to 60–75%, final analysis timeline), index.md.
