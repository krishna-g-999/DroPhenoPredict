# EmoryDrosophilaTau — what to download for cross-lab Tau replication

Verified structure (Synapse syn7274101, Emory lab — independent from BCM-DMAS/Baylor). Same DUA
pattern as BCM-DMAS: open the page, accept terms if prompted, download via your browser.

## Recommended (directly usable, small)

| File | Synapse page | Why |
|---|---|---|
| `EmoryDrosophilaTau_biospecimen_metadata.csv` | https://www.synapse.org/Synapse:syn18879642 | sample→genotype/condition key (needed to interpret everything else) |
| `EmoryDrosophilaTau_Proteomics_metadata.csv` | https://www.synapse.org/Synapse:syn18823203 | sample→channel/condition key for the proteomics |
| `baylor_proteomics_fly_proteinoutput.txt` | https://www.synapse.org/Synapse:syn7421713 | **the actual protein data** — same kind of processed table we already built a full analysis pipeline for (Findings 17). Directly comparable, real cross-lab Tau proteomics. |

Put all three in `D:\DroPhenoPredict\data\` (same as before) — I'll take it from there.

## Optional, bigger undertaking (only if you want cell-type resolution)

**scRNA-seq** — standard 10x Genomics format (`barcodes.tsv.gz` + `features.tsv.gz` + `matrix.mtx.gz`
per sample, multiple samples e.g. `UU20893_*`, `UU20894_*`, ...). This is genuinely usable but is a
bigger analysis (needs `scanpy`/cell-type annotation, not a quick add-on). List the sample files at
https://www.synapse.org/Synapse:syn41894625 if you want to pursue this — I'd recommend deciding after
we see what the proteomics/RNA cross-lab comparison shows first.

## Not currently actionable — do not download

**Bulk RNA-seq** is deposited as **raw FASTQ/BAM only** (no processed counts matrix like BCM-DMAS
had). Using it would require running our own alignment + quantification pipeline — a substantial new
piece of infrastructure, out of scope unless you specifically want to invest in it. Skip for now.
