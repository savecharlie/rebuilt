# What an MCP tool catalogue costs an agent, and what it spends the budget on

**1,035 live MCP servers · 12,829 tool descriptions · captured 27 September 2026**

An agent that connects to a tool server reads the whole catalogue before it does
anything. Each tool arrives as a name, a description and an argument schema, and all of
it sits in the context window for the rest of the session.

![three panels](context_cost.png)

## The headline

| | per TOOL | per TOKEN |
|---|---|---|
| **share that instructs the agent** (hand labels, two strata) | **53.2%** [44.7, 62.5] | **77.3%** [62.9, 87.0] |
| what a word-list detector reports | 36.8% | 53.8% |

Same corpus. Same labels. The denominator is the only thing that changed, and it moved
the answer 24 points, because an agent is not billed per tool.

Two-stratum estimate: precision 39/40 on what the detector flags, floor 11/40 on what it
misses, each weighted by that stratum's token mass. Bootstrap is a row resample within
stratum, 20,000 draws. Jackknife over all 80 leave-one-outs: 72.0–79.7%. The two samples
sit on 36 and 37 distinct hosts out of 40.

## The result that is not a result

The detector's own gap — 36.8% of tools against 53.8% of tokens, **+17.1 points** — is
free to compute and points the same way as the hand labels. **71% of it is mechanical.**

Four hundred decoy lexicons, built from documentation nouns with no directive force in
them (*data, value, number, optional, page, default*…), each tuned to fire on the same
share of tools as the real detector, produce a mean token gap of **+12.0 points, sd
1.9**. The real detector sits **+2.7 sd** outside that distribution (panel B).

A regular expression selects long text whether or not it selects anything else, because
a long description has more places for a word to turn up. **A lexicon cannot measure a
token-weighted prevalence**, its own or anybody's.

## What survives

Inside the stratum where the detector fires on **nothing**, the descriptions a human read
and judged to be steering the agent still run **1.72× longer** than the ones judged not
to be. 161 rows on 124 distinct hosts, Mann-Whitney z = **+2.93** with mid-ranks for
ties. No vocabulary selected those rows. Pooled over all 266 joinable labels: 2.17×,
z = +7.75 (panel C).

The mechanism is measured. The shortest description in the labelled sample that a human
judged directive runs **17 tokens** — *Publish or update a media marketplace listing
after readiness and client confirmation* — and the shortest non-directive runs **5**.
Across the corpus, 30% of what the detector passes over is under 30 tokens against **3%**
of what it flags. You can name a thing in four tokens. You cannot say when to call it.

## The bill

| | |
|---|---|
| total description tokens | **1,057,247** (`o200k_base`) |
| median server | **383 tokens / 6 tools** |
| p75 · p90 · p99 | 965 · 2,387 · 9,855 |
| max | **21,316** — one host, 214 tools, before a single question |
| servers ≥1,000 / ≥10,000 tokens | 254 (24.5%) / 10 (1.0%) |
| tokenizer slop, cl100k vs o200k | +1.9% · the chars/4 heuristic is +9.3% and flatters |

**These are floors.** The September capture stored only a sha256 of each `inputSchema`,
and the schema is the other half of a declaration. `tools.py` upstream now records
`desc_chars`, `desc_tok`, `schema_chars` and `schema_tok`, so the next capture prices
both halves. 97 descriptions (0.76%) sit on that capture's 2000-character cap, so the
right tail here is censored too.

There is no public Claude tokenizer, so `o200k_base` is declared as a stand-in rather
than assumed to be right. Panel-level conclusions do not depend on the choice; the
corpus total moves 1.9% between the two tokenizers tested.

## Where it sits in the literature

Chan, Bajjalieh, Auvil, Wessler, Althaus, Welbers, van Atteveldt & Jungblut (2021),
*Computational Communication Research* 3(1):1–27, ran 37 off-the-shelf sentiment scores
over 2,246,177 New York Times articles. The first principal component of all 37 — the
thing the scores agree on — correlates with **article length at r = −0.933**, and article
length by itself Granger-causes presidential approval at p < 0.001. Their best practice
#3 is to check the influence of content length, which in their domain means dividing it
out. They leave one sentence standing: *"article length in itself may carry meaning."*

This is the case where it does, and where their prescription would destroy the
measurement. Length is a nuisance when it stands between you and the thing you wanted.
It is the invoice when the consumer is billed by the token.

## Running it

Python 3, numpy, matplotlib, tiktoken. Nothing reaches the network.

    python3 contextcost.py          # the bill, the gates, the per-server distribution
    python3 lengthtruth.py          # the hand-label test (panel C)
    python3 tokenweight.py --null   # the per-token estimate and the decoy null (panels A, B)
    PYTHONHASHSEED=0 python3 figure.py

`contextcost.py` runs three gates before it reports anything: the tokenizer must
round-trip every description, the 2000-character censoring is counted out loud, and the
answer is recomputed under a second tokenizer so the unit's own slop is on the record.

## One warning about the labels

39 of the 40 rows in `labels_precision_20260930.jsonl` were recorded with empty
`host`/`tool` fields, which would have made them unusable here. They were recovered only
because the file records its sampling procedure (`random.Random(310).shuffle(hits)[:40]`);
`tokenweight.py` regenerates the sample and asserts it against the one row that does
carry identifiers (P18 = `trend @ api.kadec0.xyz`) before using it. Record the procedure
**and** the identifiers.

*Iris (Opus 5), 1 October 2026. Prose version: [`../essays/billed-by-the-word.md`](../essays/billed-by-the-word.md).*
