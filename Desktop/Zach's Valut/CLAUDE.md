# REGAL Trial Wiki — Schema & Operating Instructions

## Purpose
This wiki supports deep analysis of the REGAL Phase 3 trial (galinpepimut-S vs. best available therapy in AML CR2) with the goal of building a well-evidenced prediction of the trial's likely outcome. It is maintained entirely by the LLM; the human curates sources and directs analysis.

## Directory Layout

```
raw/
  papers/      ← immutable PDFs of clinical papers, trial documents
  reddit/      ← reddit post text, community analyses (saved as .md or .txt)
  other/       ← chat exports, news items, misc
wiki/
  index.md     ← master catalog of all wiki pages (LLM updates on every ingest)
  log.md       ← append-only activity log (LLM appends on every operation)
  overview.md  ← evolving synthesis & outcome prediction thesis
  sources/     ← one page per raw source
  concepts/    ← statistical methods, trial design concepts, endpoints
  entities/    ← GPS, REGAL trial, AML CR2, BAT arm, Sellas Life Sciences, etc.
  analyses/    ← specific quantitative analyses, scenario models, predictions
```

**Rule:** Never modify files under `raw/`. Read them, never write them.

## Page Formats

### Source pages (`wiki/sources/SLUG.md`)
```markdown
---
type: source
title: "Full title of source"
source_type: paper | reddit | chat | other
date_ingested: YYYY-MM-DD
original_file: raw/papers/filename.pdf
tags: [aml, gps, regal, statistics, ...]
---

## Summary
2-4 paragraph synthesis of key content.

## Key Claims
- Bullet list of specific factual claims with page/section refs where possible.

## Relevance to REGAL Outcome Prediction
How does this source update our probability estimate or model?

## Contradictions / Tensions
Any conflicts with other wiki pages? Name the pages.

## Open Questions
What does this source leave unresolved?
```

### Entity pages (`wiki/entities/SLUG.md`)
```markdown
---
type: entity
name: "Entity name"
aliases: [list of alternate names]
---

## Overview
What this entity is.

## Key Facts
Bulleted facts with source citations [[source-slug]].

## Relevance to REGAL
Why this entity matters for the trial prediction.

## Open Questions
```

### Concept pages (`wiki/concepts/SLUG.md`)
```markdown
---
type: concept
name: "Concept name"
---

## Definition
## How It Applies to REGAL
## Key Sources
## Notes / Caveats
```

### Analysis pages (`wiki/analyses/SLUG.md`)
```markdown
---
type: analysis
title: "Analysis title"
date: YYYY-MM-DD
sources_used: [list of source slugs]
---

## Question
## Method
## Findings
## Confidence Level
low | medium | high — with justification
## Implications for Outcome Prediction
## Caveats
```

## Key Entities to Maintain
Always keep these pages current as new sources arrive:
- `entities/gps.md` — Galinpepimut-S (mechanism, trial history, safety)
- `entities/regal-trial.md` — REGAL trial design, arms, endpoints, timeline
- `entities/aml-cr2.md` — AML second complete remission, patient population
- `entities/bat-arm.md` — Best Available Therapy comparator details
- `entities/sellas.md` — Sellas Life Sciences corporate context

## Workflows

### Ingest a new source
1. Read the source file completely.
2. Discuss key takeaways with the user before writing.
3. Write a source page in `wiki/sources/`.
4. Update or create any entity pages touched by the source.
5. Update or create any concept pages touched by the source.
6. Update `wiki/overview.md` if the source materially changes the prediction picture.
7. Update `wiki/index.md` — add entries for all new/updated pages.
8. Append to `wiki/log.md`: `## [YYYY-MM-DD] ingest | Source Title`

### Answer a query
1. Read `wiki/index.md` to find relevant pages.
2. Read those pages. Read additional pages as needed.
3. Synthesize an answer with citations to wiki pages (e.g. [[regal-trial]]).
4. If the answer is substantive, offer to file it as a new analysis page.
5. Append to `wiki/log.md`: `## [YYYY-MM-DD] query | Brief question summary`

### Lint the wiki
1. Scan all wiki pages for contradictions between pages.
2. Find orphan pages (no inbound links from other wiki pages).
3. Find concepts or entities mentioned but lacking their own page.
4. Check `overview.md` — does it still reflect the current evidence?
5. Suggest new sources to find (papers, analyses) that would reduce uncertainty.
6. Append to `wiki/log.md`: `## [YYYY-MM-DD] lint | Summary of issues found`

## Citation Convention
- Within wiki pages, link to other wiki pages as `[[slug]]` (Obsidian wikilinks).
- Cite raw sources as `[[sources/slug]]`.
- When making probabilistic claims, always state the confidence level and what evidence would change it.

## Overview Page Conventions
`overview.md` is the most important page. It should always contain:
1. **Current prediction**: a stated probability range for GPS meeting primary endpoint, with reasoning.
2. **Key evidence for**: bulleted, cited.
3. **Key evidence against**: bulleted, cited.
4. **Critical uncertainties**: what we don't know that matters most.
5. **What would change the prediction**: specific evidence that would shift the estimate significantly.

## Log Format
Each log entry starts with `## [YYYY-MM-DD] operation | description` so entries are grep-parseable.
