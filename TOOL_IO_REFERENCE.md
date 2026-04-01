# Primer Design Pipeline — Tool Input/Output Reference

> **Pipeline:** Progressive refinement — 20,000 raw → filtered → checked → final assay sets  
> **Organism:** SARS-CoV-2 (5 genomes: MN908947.3, MN985325.1, MN988713.1, MN938384.1, MN975262.1)

---

## Pipeline Flow Overview

```
test5.fasta (5 genomes, 29,903 bp each)
    │
    ├─ ALIGN ─────────────────── MAFFT ──────────────────► mafft_aligned.fasta (MSA)
    │                            Clustal Omega ──────────► clustalo_aligned.fasta (MSA)
    │
    ├─ GENERATE & PRE-FILTER ─── primer3_core ───────────► 10,000 LEFT + 10,000 RIGHT = 20,000 raw primers
    │                            primer3-py filter ──────► 1,061 filtered primers (TSV + FASTA + pairs)
    │                            primer3-py (multi) ─────► per-genome primer pairs
    │
    ├─ FILTER & SCORE ────────── MELTING ────────────────► nearest-neighbour Tm (50 sample)
    │                            oligo-melting ──────────► Tm all filtered primers
    │                            BLAST ──────────────────► off-target hit table
    │                            MFEprimer ──────────────► specificity + dimer + hairpin reports
    │                            seqkit ─────────────────► FASTA stats
    │                            bowtie2 + samtools ─────► alignment SAM
    │                            isPcr ──────────────────► amplicons FASTA (exact PCR)
    │                            ipcress ────────────────► amplimers (mismatch-tolerant)
    │                            primersearch ───────────► EMBOSS amplimers (mismatch-tolerant)
    │                            tntblast ───────────────► thermodynamic PCR simulation
    │
    ├─ OPTIMIZE / PANEL SELECTION ─ olivar ──────────────► SADDLE-optimised tiling (3 kb demo)
    │                               NGS-PrimerPlex ──────► multiplex panel (2 target regions)
    │
    └─ SPECIALIZED WORKFLOWS ──── PrimalScheme3 ─────────► tiling scheme BED + HTML
                                 varvamp ───────────────► variation-aware tiling scheme
                                 PUPpy ─────────────────► taxon-specific primers (E. coli vs Salmonella)
                                 DegePrime ─────────────► degenerate primer windows from MSA
                                 PrimerServer2 ─────────► target-focused primer pairs (5 regions)
                                 primerdiffer ──────────► discriminatory primers (10 genome pairs)
```

---
## ALIGN

### 1a) MAFFT

**Purpose:** Multiple sequence alignment of all 5 genomes — required by DegePrime, varvamp, PrimalScheme3.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `test5.fasta` | FASTA | 5 unaligned SARS-CoV-2 genomes (~29,900 bp each) |
| **OUT** | `mafft_aligned.fasta` | FASTA | 5 aligned genomes (same length, gaps inserted) |

**Command:** `mafft --auto --thread 4 test5.fasta`

**Output sample (`mafft_aligned.fasta`):**
```
>MN908947.3 Severe acute respiratory syndrome coronavirus 2 ...
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT
...
acagtgaacaatgctagggagagctgcctatatggaagagccctaatgtgtaaaattaat  ← last genome row
aaaaaaaaaaa------------                                        ← trailing gaps from alignment
```

---

---

### 1b) Clustal Omega

**Purpose:** Alternative MSA — for comparison/validation.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `test5.fasta` | FASTA | 5 unaligned genomes |
| **OUT** | `clustalo_aligned.fasta` | FASTA | 5 aligned genomes |

**Command:** `clustalo -i test5.fasta -o clustalo_aligned.fasta --outfmt=fasta --force`

**Output sample (`clustalo_aligned.fasta`):**
```
>MN908947.3 Severe acute respiratory syndrome coronavirus 2 ...
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT
...
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa--  ← last line, gap-padded
```

---

## STAGE 2 — Raw Candidate Generation

---

## GENERATE & PRE-FILTER

### 2a) primer3_core

**Purpose:** Generate maximum possible primer candidates — cast wide net, filter later.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `primer3_input_left.txt` | primer3 boulder-IO | MN908947.3 full genome (29,903 bp), `PRIMER_NUM_RETURN=10000`, `PRIMER_TASK=pick_primer_list` |
| **IN** | `primer3_input_right.txt` | primer3 boulder-IO | Same genome, right primers only |
| **OUT** | `primer3_output_left.txt` | primer3 boulder-IO | 10,000 left primer sequences with Tm, GC, position |
| **OUT** | `primer3_output_right.txt` | primer3 boulder-IO | 10,000 right primer sequences |
| **OUT** | `primer3_output.txt` | primer3 boulder-IO | Combined (left + right concatenated) |

**Input parameters:**
```
PRIMER_TASK=pick_primer_list
PRIMER_OPT_SIZE=20  PRIMER_MIN_SIZE=18  PRIMER_MAX_SIZE=25
PRIMER_OPT_TM=60.0  PRIMER_MIN_TM=57.0  PRIMER_MAX_TM=63.0
PRIMER_MIN_GC=40.0  PRIMER_MAX_GC=60.0
PRIMER_NUM_RETURN=10000
PRIMER_PRODUCT_SIZE_RANGE=100-400
```

**Output sample (boulder-IO key=value, one record per primer):**
```
SEQUENCE_ID=MN908947.3_LEFT
PRIMER_LEFT_0_PENALTY=0.028505
PRIMER_LEFT_0_SEQUENCE=GCCGCTGTTGATGCACTATG
PRIMER_LEFT_0=17169,20
PRIMER_LEFT_0_TM=59.971
...
PRIMER_LEFT_9999_SEQUENCE=...  ← up to 10,000 entries
=
```
Right file equivalent:
```
PRIMER_RIGHT_0_SEQUENCE=CATAGTGCATCAACAGCGGC
PRIMER_RIGHT_0=17188,20
PRIMER_RIGHT_0_TM=59.971
```

**Result:** 10,000 LEFT + 10,000 RIGHT = **20,000 raw primers**

---

## STAGE 3a — Pre-filter: Tm + Structure

---

### 3a) primer3-py (filter pass)

**Purpose:** Fast cheap filters — reject bad Tm, hairpins, homodimers from all 20,000 raw primers.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `primer3_output_left.txt` | primer3 boulder-IO | 10,000 left primers |
| **IN** | `primer3_output_right.txt` | primer3 boulder-IO | 10,000 right primers |
| **OUT** | `all_primers_tm.tsv` | TSV | Tm + structure for ALL ~20,000 unique primers |
| **OUT** | `primers_filtered.tsv` | TSV | Primers passing all filters |
| **OUT** | `primers_filtered.fasta` | FASTA | Same primers in FASTA format |
| **OUT** | `pairs_filtered.tsv` | TSV | Matched F+R pairs from filtered set |
| **OUT** | `pairs_100_ispcr.txt` | tab-delimited | Top 100 pairs for isPcr |
| **OUT** | `pairs_100_ipcress.txt` | space-delimited | Top 100 pairs for ipcress |
| **OUT** | `pairs_100_primersearch.txt` | tab-delimited | Top 100 pairs for primersearch |
| **OUT** | `pairs_100_tntblast.fasta` | FASTA | Top 100 pairs for tntblast |
| **OUT** | `primers_melting_50.txt` | plain text | 50 primer sequences for MELTING |

**Filters applied:** Tm 57–63°C AND hairpin dG > −2000 cal/mol AND homodimer dG > −2000 cal/mol

**`all_primers_tm.tsv` sample:**
```
primer_seq              Tm_C    hairpin_dG  housedimer_dG
TCAGTTACGTGCCAGATCAG    60.19   0           -4723
GCAACTGAGGGAGCCTTGA     62.65   253         -2419
AGGCAGGTCCTTGATGTCACA   64.25   -381        -2955   ← FAILS Tm > 63
```

**`primers_filtered.tsv` sample:**
```
primer_seq                  Tm_C    hairpin_dG  housedimer_dG
ACAACAGCCCTTGAGACAACTA      62.30   0           -1914
TTACCAACCACCACAAACCT        60.12   0           651
CAAGCCTCTTCTCGTTCCTCAT      62.79   0           -1902
```

**`primers_filtered.fasta` sample:**
```
>primer_0
ACAACAGCCCTTGAGACAACTA
>primer_1
TTACCAACCACCACAAACCT
```

**seqkit stats on filtered FASTA:**
```
num_seqs=1,061  sum_len=22,063  min_len=18  avg_len=20.8  max_len=23
```

**`pairs_100_ispcr.txt` format:** `name<TAB>F_seq<TAB>R_seq`
```
pair1   GAAGTGGGTTTTGTCGTGCC    CGCCAACAATAAGCCATCCG
pair2   GGTTCACCGCTCTCACTCAA    TAGCCCATCTGCCTTGTGTG
```

**`pairs_100_ipcress.txt` format:** `name F_seq R_seq min_size max_size`
```
pair1 GAAGTGGGTTTTGTCGTGCC CGCCAACAATAAGCCATCCG 50 500
```

**Result:** 20,000 → **1,061 primers** pass; matched into pairs for downstream tools

---

## STAGE 3b — Thermodynamics on Filtered Primers

---

### 7c) primer3-py (multi-genome)

**Purpose:** Direct Python API primer design from all 5 genomes independently.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `test5.fasta` | FASTA | All 5 SARS-CoV-2 genomes |
| **OUT** | `primer3py_output.txt` | TSV | genome_id, pairs_returned, cleaned_ambiguous_bases |

**Parameters:** size 18–25 bp, Tm 57–63°C, GC 40–60%, 50 pairs returned per genome, product 100–400 bp  
**Note:** Ambiguous bases (IUPAC non-ACGTN) replaced with N before design.

**Output sample (`primer3py_output.txt`):**
```
genome_id    pairs_returned  cleaned_ambiguous_bases
MN908947.3   50              0
MN985325.1   50              0
MN988713.1   50              8
MN938384.1   50              0
MN975262.1   50              0   ← last genome
```

---

## Summary Table

| Tool | Stage | Input | Key Output | Role |
|------|-------|-------|------------|------|
| MAFFT | 1 | `test5.fasta` | `mafft_aligned.fasta` | MSA |
| Clustal Omega | 1 | `test5.fasta` | `clustalo_aligned.fasta` | MSA (alt) |
| primer3_core | 2 | genome + params | `primer3_output_{left,right}.txt` | 20k raw primers |
| primer3-py | 3a | 20k raw primers | `primers_filtered.{tsv,fasta}` | Tm+structure filter |
| MELTING | 3b | 50 primer seqs | `melting_results.txt` | ΔH/ΔS/Tm (NN model) |
| oligo-melting | 3b | 1,061 filtered primers | `oligomelting_results.txt` | Tm at 1M + 50mM Na⁺ |
| blastn | 4a | filtered primers vs 5 genomes | `blast_results.txt` | Off-target hits |
| MFEprimer | 4a | filtered primers | `mfeprimer_{spec,dimer,hairpin}.txt` | Specificity + structure |
| primerdiffer | 4a | genome pairs | `primerdiffer_all/` | Discriminatory primers |
| seqkit | 4a | filtered FASTA | `seqkit_primer_stats.txt` | FASTA stats |
| bowtie2 | 4a | filtered primers vs 5 genomes | `bowtie2_primers.sam` | Alignment specificity |
| isPcr | 4b | top 100 pairs | `ispcr_output.fasta` | Exact amplicons |
| ipcress | 4b | top 100 pairs | `ipcress_output.txt` | Mismatch-tolerant PCR |
| primersearch | 4b | top 100 pairs | `primersearch_output.txt` | EMBOSS PCR (10% mm) |
| tntblast | 4b | 1 pair | `tntblast_output.txt` | Thermodynamic PCR sim |
| PrimalScheme3 | 5 | MSA | `ps3_out/primer.bed` + `amplicon.bed` | Tiling scheme |
| varvamp | 5 | MSA | `varvamp_out/primers.tsv` + BED | Variation-aware tiling |
| olivar | 5 | 3kb FASTA | `olivar_out/sars2_tiling.csv` + BED | SADDLE tiling |
| PUPpy | 6 | CDS target + nontarget dirs | `ResultDB.tsv` → `UniquePrimerTable.tsv` | Taxon-specific |
| NGS-PrimerPlex | 6 | BED targets + genome | `*_combination_1.fa` | Multiplex panel |
| DegePrime | 7 | MSA | `degeprime_output.tsv` | Degenerate primers |
| PrimerServer2 | 7 | bracket-format FASTA | TSV primer pairs per region | Target-focused design |
| primer3-py (multi) | 7 | `test5.fasta` | `primer3py_output.txt` | Per-genome design |

---

## FILTER & SCORE

### 3b-i) MELTING (Java)

**Purpose:** Nearest-neighbour thermodynamic Tm with Na⁺ salt correction.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `primers_melting_50.txt` | plain text | 50 primer sequences (one per line) |
| **OUT** | `melting_out/melting_results.txt` | text | ΔH, ΔS, Tm for each primer |

**Command per primer:** `java -jar melting5.jar -S <SEQ> -H dnadna -P 0.00000025 -E Na=0.05`

**Output sample:**
```
The MELTING results are :
Enthalpy : -150,000 cal/mol ( -627,000 J/mol)
Entropy  : -405.3 cal/mol-K ( -1,694.15 J/mol-K)
Melting temperature : 52.91 degrees C.
```

**Note:** Nearest-neighbour model, [Na⁺]=50 mM, [oligo]=250 nM

---

---

### 3b-ii) oligo-melting

**Purpose:** Python-based Tm calculation on ALL filtered primers with salt correction.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `primers_filtered.fasta` | FASTA | 1,061 filtered primers |
| **OUT** | `oligomelting_results.txt` | text | Tm at 1M Na⁺, Tm at 50mM Na⁺, ΔH, ΔS per primer |

**Command per primer:** `oligomelting <SEQ> --conc-na 50`

**Output sample:**
```
Sequence:    AGGTTTGTGGTGGTTGGTAA
Tm (1M Na+): 54.91 C
Tm (50mM Na+): 40.03 C
dH: -150.0 cal/mol
dS: -0.405300 cal/mol/K
```

---

## STAGE 4a — Specificity on Filtered Primers

---

### 4a-i) makeblastdb + blastn

**Purpose:** Check off-target binding of filtered primers across all 5 SARS-CoV-2 genomes.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `primers_filtered.fasta` | FASTA | 1,061 filtered primers as query |
| **IN** | `test5.fasta` → `sars2_blastdb.*` | BLAST DB | 5-genome BLAST database |
| **OUT** | `blast_results.txt` | tabular (fmt 6) | Hit table: qseqid, sseqid, pident, length, mismatches, gapopen, qstart, qend, sstart, send, evalue, bitscore |
| **LOG** | `makeblastdb.log` | text | DB creation stats |

**Command:**
```bash
makeblastdb -in test5.fasta -dbtype nucl -out sars2_blastdb
blastn -query primers_filtered.fasta -db sars2_blastdb -task blastn-short \
       -word_size 7 -outfmt 6 -out blast_results.txt
```

**Output sample:**
```
primer_0  MN975262.1  100.000  22  0  0  1   22  25306  25285  9.47e-08  44.1
primer_0  MN975262.1  100.000   9  0  0  5   13  5378   5370   5.4       18.3
primer_0  MN975262.1  100.000   9  0  0  14  22  17353  17345  5.4       18.3
```

---

---

### 4a-ii) MFEprimer

**Purpose:** Binding specificity, primer dimers, and hairpin thermodynamics on filtered primers.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `primers_filtered.fasta` | FASTA | 1,061 filtered primers |
| **IN** | `test5.fasta` (indexed) | MFEprimer index | 5-genome reference |
| **OUT** | `mfeprimer_spec.txt` | text report | Per-primer: Length, GC%, Tm, ΔG, binding counts (plus/minus strand) |
| **OUT** | `mfeprimer_spec.txt.spec.tsv` | TSV | Amplicon table: amp coordinates, primer sequences, Tm, GC, ΔG, genome name |
| **OUT** | `mfeprimer_dimer.txt` | text report | Dimer analysis: primer sequence, Length, GC%, Tm, ΔG |
| **OUT** | `mfeprimer_hairpin.txt` | text report | Hairpin analysis: same columns as dimer |

**Commands:**
```bash
mfeprimer index -i test5.fasta
mfeprimer spec    -i primers_filtered.fasta -d test5.fasta -o mfeprimer_spec.txt
mfeprimer dimer   -i primers_filtered.fasta -o mfeprimer_dimer.txt
mfeprimer hairpin -i primers_filtered.fasta -o mfeprimer_hairpin.txt
```

**`mfeprimer_spec.txt` sample:**
```
Primer ID    Sequence (5'→3')              Length  GC%    Tm(°C)  Dg(kcal/mol)  +Bind  -Bind
primer_0     AGGTTTGTGGTGGTTGGTAA          20      45.00  57.67   -20.51         0      10
primer_1     TGATGGCTACCCTCTTGAGT          20      50.00  58.27   -20.80         5       0
primer_8     CAAGGCGTTCCAATTAACACCA        22      45.45  60.22   -23.03        10       0
```

**`mfeprimer_spec.txt.spec.tsv` sample:**
```
Amp_1  MN985325.1  28481  29561  46.81  1081  primer_8  CAAGGCGTTCCAATTAACACCA  ...  primer_821  AGACCACACAAGGCAGATGG
```

**`mfeprimer_dimer.txt` / `mfeprimer_hairpin.txt` sample:**
```
primer_0  ACAACAGCCCTTGAGACAACTA  22  45.45  60.09  -22.75
primer_3  CAAGCCTCTTCTCGTTCCTCAT  22  50.00  60.62  -23.20
```

---

---

### 4a-iv) seqkit

**Purpose:** Quick statistics on the filtered primer FASTA.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `primers_filtered.fasta` | FASTA | 1,061 filtered primers |
| **OUT** | `seqkit_primer_stats.txt` | TSV | file, format, type, num_seqs, sum_len, min_len, avg_len, max_len |

**Output:**
```
file                               format  type  num_seqs  sum_len  min_len  avg_len  max_len
primers_filtered.fasta             FASTA   DNA    1,061    22,063      18     20.8       23
```

---

---

### 4a-v) bowtie2 + samtools

**Purpose:** Alignment-based specificity — maps all filtered primers against all 5 genomes.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `test5.fasta` → `sars2_bt2.*` | bowtie2 index | 5-genome index (.bt2 files) |
| **IN** | `primers_filtered.fasta` | FASTA | 1,061 primers, treated as single-end reads |
| **OUT** | `bowtie2_primers.sam` | SAM | Alignment hits (header + aligned reads) |
| **LOG** | `bowtie2_build.log` | text | Index build settings |
| **LOG** | `bowtie2.log` | text | Alignment summary (% aligned) |

**Command:**
```bash
bowtie2-build test5.fasta sars2_bt2
bowtie2 -x sars2_bt2 -f -U primers_filtered.fasta --no-unal -N 1 -L 18 -S bowtie2_primers.sam
```

**SAM header sample:**
```
@SQ  SN:MN908947.3  LN:29903
@SQ  SN:MN985325.1  LN:29882
...
primer_0  16  MN938384.1  25253  1  22M  *  0  0  TAGTTGTCTCAAGGGCTGTTGT  ...  AS:i:0  NM:i:0
```

---

## STAGE 4b — In-silico PCR on Top 100 Pairs

---

### 4b-i) isPcr

**Purpose:** Exact in-silico PCR simulation — finds amplicons from exact primer matches.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `test5.fasta` | FASTA | All 5 genomes as template |
| **IN** | `pairs_100_ispcr.txt` | tab-delimited | `name<TAB>F_seq<TAB>R_seq` — 100 pairs |
| **OUT** | `ispcr_output.fasta` | FASTA | Amplicon sequences with coordinates in header |

**Command:** `isPcr test5.fasta pairs_100_ispcr.txt stdout`

**Output sample:**
```
>MN985325.1:28428+29565 pair2 1138bp GGTTCACCGCTCTCACTCAA TAGCCCATCTGCCTTGTGTG
GGTTCACCGCTCTCACTCAA...catggcaaggaagaccttaaa...TAGCCCATCTGCCTTGTGTG
```

Header format: `genome:start+end pairname size F_seq R_seq`

---

---

### 4b-ii) ipcress

**Purpose:** Mismatch-tolerant in-silico PCR (up to 2 mismatches).

| | File | Format | Content |
|---|---|---|---|
| **IN** | `pairs_100_ipcress.txt` | space-delimited | `name F_seq R_seq min_size max_size` |
| **IN** | `test5.fasta` | FASTA | Template genomes |
| **OUT** | `ipcress_output.txt` | text | Amplimers found; "-- completed ipcress analysis" if no hits |

**Command:** `ipcress --input pairs_100_ipcress.txt --sequence test5.fasta --mismatch 2`

**Output (this run — no amplimers at 2mm threshold):**
```
-- completed ipcress analysis
```
**Expected format when amplimers found:**
```
ipcress: MN908947.3 pair1 350 A 8846 0 B 9196 0 forward
```
Fields: `ipcress: seqname experiment_name product_len strand fwd_pos fwd_mm rev_pos rev_mm orientation`

---

---

### 4b-iii) primersearch (EMBOSS)

**Purpose:** EMBOSS mismatch-tolerant PCR simulation (10% mismatch allowed).

| | File | Format | Content |
|---|---|---|---|
| **IN** | `pairs_100_primersearch.txt` | tab-delimited | `name<TAB>F_seq<TAB>R_seq` |
| **IN** | `test5.fasta` | FASTA | Template |
| **OUT** | `primersearch_output.txt` | text | Amplimer blocks with hit positions and mismatch counts |

**Command:**
```bash
primersearch -seqall test5.fasta -infile pairs_100_primersearch.txt \
             -mismatchpercent 10 -outfile primersearch_output.txt
```

**Output sample:**
```
Primer name pair1
Amplimer 1
    Sequence: MN908947.3
    GAAGTGGGTTTTGTCGTGCC hits forward strand at 8846 with 0 mismatches
    CGCCAACAATAAGCCATCCG hits reverse strand at [4364] with 0 mismatches
    Amplimer length: 16695 bp
```

---

---

### 4b-iv) tntblast

**Purpose:** Full thermodynamic PCR simulation — reports Tm, ΔG, hairpin, homodimer, heterodimer.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `tntblast_assay.txt` | TSV | `#name<TAB>forward<TAB>reverse` — 1 hardcoded pair |
| **IN** | `test5.fasta` | FASTA | Template |
| **OUT** | `tntblast_output.txt` | text | Full thermodynamic report per pair |
| **LOG** | `tntblast.log` | text | Run output |

**Assay file:**
```
#name   forward                 reverse
pair1   TCGAACTGCACCTCATGGTC    GACTTTAGATCGGCGCCGTA
```

**Output sample:**
```
name = pair1
forward primer = 5' TCGAACTGCACCTCATGGTC 3'
reverse primer = 5' GACTTTAGATCGGCGCCGTA 3'
forward primer tm = 59.8755
reverse primer tm = 59.9081
forward primer hairpin tm = 25.1165
reverse primer hairpin tm = 45.9363
forward primer homodimer tm = 0
heterodimer tm = 0
forward primer dG[-19.3436] = dH[-156.7] - T*dS[-0.442871]
forward primer mismatches = 0
forward primer %GC = 55
```

---

## STAGE 5 — Tiling Schemes

---

## OPTIMIZE / PANEL SELECTION

### 5c) olivar

**Purpose:** SADDLE-optimization tiling scheme — risk-minimised primers (first 3,000 bp demo).

| | File | Format | Content |
|---|---|---|---|
| **IN** | `olivar_short.fasta` | FASTA | First 3,000 bp of MN908947.3 |
| **OUT** | `olivar_db/sars2_olivar.olvr` | binary | olivar database file |
| **OUT** | `olivar_out/sars2_tiling.csv` | CSV | amplicon_id, pool, fP, rP, start, insert_start, insert_end, end, amplicon_seq, insert_seq |
| **OUT** | `olivar_out/sars2_tiling.scheme.bed` | BED | Primer positions with pool assignment |
| **OUT** | `olivar_out/sars2_tiling.fasta` | FASTA | Primer sequences |
| **OUT** | `olivar_out/sars2_tiling.json` | JSON | Scheme metadata |
| **OUT** | `olivar_out/sars2_tiling.html` | HTML | Interactive scheme visualisation |
| **OUT** | `olivar_out/sars2_tiling_Loss.html` | HTML | SADDLE loss curve plot |
| **OUT** | `olivar_out/sars2_tiling_risk.csv` | CSV | Per-position risk scores |

**Commands:**
```bash
olivar build olivar_short.fasta --output olivar_db --title sars2_olivar --threads 4
olivar tiling olivar_db/sars2_olivar.olvr --output olivar_out --title sars2_tiling \
              --max-amp-len 400 --seed 42 --threads 4
```

**`sars2_tiling.csv` sample:**
```
amplicon_id,      pool,fP,                              rP,                  start,insert_start,insert_end,end
sars2_tiling_1,   1,   ccaaccaactttcgatctcttgtag,       tcctccacggagtctccaaa, 32,   57,          355,        375
sars2_tiling_2,   2,   gtgcactcacgcagtataattaataact,    gctgttcaagttgaggcaaaacg, 114, 142,       441,        464
```

**`sars2_tiling.scheme.bed` sample:**
```
input-seq-1  32   56   sars2_tiling_1_LEFT   1  +  CCAACCAACTTTCGATCTCTTGTAG
input-seq-1  356  375  sars2_tiling_1_RIGHT  1  -  TCCTCCACGGAGTCTCCAAA
```

---

## STAGE 6 — Taxon-specific + Multiplex

---

### 6b) NGS-PrimerPlex

**Purpose:** Design multiplexed primer panel for 2 SARS-CoV-2 target regions (amplicon-based NGS).

| | File | Format | Content |
|---|---|---|---|
| **IN** | `npp_targets.bed` | BED | 2 target regions: `MN908947.3 400 900 target1`, `MN908947.3 1200 1700 target2` |
| **IN** | `test5.fasta` | FASTA | Reference genome |
| **OUT** | `npp_targets_primers_combination_1.fa` | FASTA | Final primer sequences |
| **OUT** | `npp_targets_primers_combination_1_info.xls` | XLS | Primer info table |
| **OUT** | `npp_targets_primers_combination_1_internal_amplicons.fa` | FASTA | Internal amplicon sequences |
| **OUT** | `npp_targets_all_draft_primers.xls` | XLS | All candidate primers before filtering |
| **LOG** | `npp_run.log` | text | Full run log with parameter dump |
| **LOG** | `npp_targets.log` | text | Design progress (% complete per region) |

**Key parameters:**
```
--min-amplicon-length 100  --max-amplicon-length 300  --optimal-amplicon-length 200
--min-primer-melting-temp 57  --max-primer-melting-temp 63  --optimal-primer-melting-temp 60
--max-primer-nonspecific 2  --primers-number1 2  --return-variants-number 1
```

**Output sample (`npp_targets_primers_combination_1.fa`):**
```
>target1_1_3_F
GGAGGAGGTCTTATCAGAGG
>target1_1_3_R
...
```
*(`.xls` files are binary Excel — open in spreadsheet app for primer info table)*

---

## STAGE 7 — Other Independent Generators

---

## SPECIALIZED WORKFLOWS

### 5a) PrimalScheme3

**Purpose:** MSA-aware amplicon tiling scheme across the full genome.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `mafft_aligned.fasta` | FASTA (MSA) | 5-genome alignment |
| **OUT** | `ps3_out/primer.bed` | BED v3 | Primer positions: chrom, start, end, name, pool, strand, sequence, pc=count |
| **OUT** | `ps3_out/amplicon.bed` | BED | Amplicon coordinates: chrom, start, end, name, pool |
| **OUT** | `ps3_out/config.json` | JSON | Full run parameters |
| **OUT** | `ps3_out/plot.html` | HTML | Interactive coverage plot |
| **OUT** | `ps3_out/primer.html` | HTML | Primer details |
| **LOG** | `ps3.log` | text | Progress (18,787 primer pairs evaluated) |

**Command:** `primalscheme3 scheme-create --msa mafft_aligned.fasta --output ps3_out --amplicon-size 400`

**Key config parameters:**
```json
{"amplicon_size": 400, "n_pools": 2, "primer_tm_min": 59.5, "primer_tm_max": 62.5,
 "primer_gc_min": 30, "primer_gc_max": 55, "primer_size_min": 19, "primer_size_max": 36,
 "version": "3.3.0"}
```

**`primer.bed` sample:**
```
MN908947.3  22   49   aeaddc3a_1_LEFT_1   1  +  GGTAACAAACCAACCAACTTTCGATCT  pc=4
MN908947.3  462  487  aeaddc3a_1_RIGHT_1  1  -  CGAACGTTTGATGAACACATAGGGC    pc=5
MN908947.3  421  449  aeaddc3a_2_LEFT_1   2  +  TTAGTAGAAGTTGAAAAAGGCGTTTTGC pc=5
```

**`amplicon.bed` sample:**
```
MN908947.3  22   487  aeaddc3a_1  1
MN908947.3  421  886  aeaddc3a_2  2
MN908947.3  822  1285 aeaddc3a_3  1
```

---

---

### 5b) varvamp

**Purpose:** Variation-aware tiling — designs primers minimising ambiguous positions.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `mafft_aligned.fasta` | FASTA (MSA) | 5-genome alignment |
| **OUT** | `varvamp_out/primers.tsv` | TSV | amplicon_name, amplicon_length, primer_name, pool, start, stop, seq, size, gc_best, temp_best, mean_gc, mean_temp, penalty, off_target_amplicons |
| **OUT** | `varvamp_out/amplicons.bed` | BED | Amplicon coordinates |
| **OUT** | `varvamp_out/primers.bed` | BED | Primer positions |
| **OUT** | `varvamp_out/primers_pool_1.fasta` | FASTA | Pool 1 primer sequences |
| **OUT** | `varvamp_out/primers_pool_2.fasta` | FASTA | Pool 2 primer sequences |
| **OUT** | `varvamp_out/primers_to_amplicon_assignment.tabular` | tabular | Primer→amplicon mapping |
| **OUT** | `varvamp_out/amplicon_plot.pdf` | PDF | Coverage plot |
| **LOG** | `varvamp.log` | text | Run log |

**Command:** `varvamp tiled mafft_aligned.fasta varvamp_out`

**`primers.tsv` sample:**
```
amlicon_name  len   primer_name          pool  start  stop  seq                       size  gc_best  temp_best
varVAMP_0     1335  varVAMP_0_LEFT       1     28     51    CAAACCAACCAACTTTCGATCTCT  24    41.7     60.0
varVAMP_0     1335  varVAMP_0_RIGHT      1     1341   1362  ACAACAGCATTTTGGGGTAAGT    22    40.9     58.9
varVAMP_2     1088  varVAMP_2_LEFT       1     2395   2416  CACGCACTCAAAGGGATTGTAC    22    50.0     60.0
```

---

---

### 6a) PUPpy

**Purpose:** Design primers uniquely specific to E. coli K-12 (not hitting Salmonella).

**Two-step tool:**

#### Step 1: puppy-align

| | File | Format | Content |
|---|---|---|---|
| **IN** | `bacterial_cds/target/` | dir of FASTA | E. coli K-12 CDS sequences (20 genes) |
| **IN** | `bacterial_cds/nontarget/` | dir of FASTA | Salmonella Typhimurium CDS (20 genes) |
| **OUT** | `puppy_align/ResultDB.tsv` | TSV | All-vs-all MMseqs2 alignment: query, target, qlen, tlen, alnlen, qstart, qend, tstart, tend, pident, qcov, tcov, evalue |
| **OUT** | `puppy_align/align_logfile.txt` | text | MMseqs2 alignment log |

**`ResultDB.tsv` sample:**
```
query                              target                              qlen  tlen  alnlen  pident   qcov   tcov   evalue
EcoliK12-..._NP_414542.1_1        EcoliK12-..._NP_414542.1_1         66    66    66      100.000  1.000  1.000  7.3E-32
EcoliK12-..._NP_414574.1_33       SalmonellaTM-..._NP_459072.1_67    3222  3228  3219    91.000   0.999  0.997  0.000
EcoliK12-..._NP_414607.1_65       SalmonellaTM-..._NP_459110.1_105   765   768   761     82.200   0.991  0.991  5.4E-220
```

#### Step 2: puppy-primers

| | File | Format | Content |
|---|---|---|---|
| **IN** | `bacterial_cds/target/` | dir of FASTA | E. coli CDS |
| **IN** | `puppy_align/ResultDB.tsv` | TSV | Alignment results |
| **OUT** | `puppy_primers/UniquePrimerTable.tsv` | TSV | Taxon-specific primer pairs with sequences, Tm, positions |

**`UniquePrimerTable.tsv` sample:**
```
species   gene                    ...  pair_penalty  amplicon_size  F_primer              R_primer              F_tm    R_tm    F_GC   R_GC
EcoliK12  cds_NP_416485.5_1966   ...  0.062         95             GATGTACGCGCAGAAAGCTG  TCGAACAGGGCCACTTCATC  59.97   60.04   55.0   55.0
EcoliK12  cds_NP_416485.5_1966   ...  0.079         76             GAACTTCACCAGCAACGCAG  AGCCTTACCCTGTTCGTTCG  60.04   60.04   55.0   55.0
```

---

---

### 7a) DegePrime

**Purpose:** Sliding-window degenerate primer design from MSA — maximises coverage across variants.

**Two-step tool (Perl scripts):**

#### Step 1: TrimAlignment.pl

| | File | Format | Content |
|---|---|---|---|
| **IN** | `mafft_aligned.fasta` | FASTA (MSA) | 5-genome alignment |
| **OUT** | `mafft_trimmed.fasta` | FASTA (MSA) | Alignment with low-coverage columns removed (`-min 0.9`) |

#### Step 2: DegePrime.pl

| | File | Format | Content |
|---|---|---|---|
| **IN** | `mafft_trimmed.fasta` | FASTA (MSA) | Trimmed alignment |
| **OUT** | `degeprime_output.tsv` | TSV | Sliding window results |

**`degeprime_output.tsv` columns:** Pos, NumberSpanning, UniqueMers, Entropy, PrimerDeg, PrimerSeq, NumberMatching, FractionMatching

**Sample:**
```
Pos  NumberSpanning  UniqueMers  Entropy  PrimerDeg  PrimerSeq             NumberMatching  FractionMatching
20   5               1           0        1          GTAGATCTGTTCTCTAAACG  5               1
21   5               1           0        1          TAGATCTGTTCTCTAAACGA  5               1
22   5               1           0        1          AGATCTGTTCTCTAAACGAA  5               1
```

**Parameters:** `-d 12` (max degeneracy), `-l 20` (primer length)

---

---

### 7b) PrimerServer2 (primertool)

**Purpose:** Target-focused primer design for 5 SARS-CoV-2 diagnostic regions with BLAST-based specificity checking.

| | File | Format | Content |
|---|---|---|---|
| **IN** | `ps2_sars2_targets.fasta` | FASTA | 5 target regions with flanks in bracket notation: `left[TARGET]right` |
| **IN** | `test5.fasta` (BLAST-indexed) | BLAST DB | Reference for specificity |
| **OUT** | `ps2_sars2_targets.fasta.json` | JSON | Intermediate design data |
| **OUT** | Printed TSV to stdout / `primerserver2_final_results.tsv` | TSV | All primer pairs per region |
| **LOG** | `ps2.log` | text | Design progress (5 regions × primer3 + BLAST) |

**Target regions in input file** (bracket notation):
```
>spike_NTD
<500bp flank>[positions 500-600 of MN908947.3]<500bp flank>

>spike_RBD
<500bp flank>[positions 1200-1350]<500bp flank>

>orf1ab_nsp3
<500bp flank>[positions 4500-4650]<500bp flank>

>envelope
<500bp flank>[positions 26200-26300]<500bp flank>

>nucleocapsid
<500bp flank>[positions 28200-28350]<500bp flank>
```

**Command:**
```bash
primertool design ps2_sars2_targets.fasta test5.fasta \
    --primer-num-return 10000 --primer-num-retain 10000 \
    -t primerserver2_final_results.tsv
```

**Output sample (`primerserver2_final_results.tsv`):**
```

---

### 4a-iii) primerdiffer

**Purpose:** Find primers that discriminate between genome pairs (diagnostic primer design).

| | File | Format | Content |
|---|---|---|---|
| **IN** | 5 individual genome FASTA files | FASTA | `MN908947.3.fasta`, `MN985325.1.fasta`, etc. |
| **OUT** | `primerdiffer_all/<g1>_vs_<g2>/` | text files | Discriminatory primers for each of 10 genome pairs |

**Command per pair:**
```bash
primerdesign.py -g1 genome1.fasta -g2 genome2.fasta -pos g1:1-29903 -d outdir/
```

**Output sample (per pair dir, e.g. `MN908947.3_vs_MN938384.1/`):**
```
>primer_g1_discriminatory_1
GTAGATCTGTTCTCTAAACG
...
>primer_g1_discriminatory_N   ← last discriminatory primer for this pair
TTTGTCACGCACTTTCCTGT
```

**Note:** Runs all 10 combinations (C(5,2)) of 5 genomes

---

---

---

## Summary Table

| Tool | Stage | Input | Key Output | Role |
|------|-------|-------|------------|------|
| MAFFT | Align | `test5.fasta` | `mafft_aligned.fasta` | MSA |
| Clustal Omega | Align | `test5.fasta` | `clustalo_aligned.fasta` | MSA (alt) |
| primer3_core | Generate & Pre-filter | genome + params | `primer3_output_{left,right}.txt` | 20k raw primers |
| primer3-py | Generate & Pre-filter | 20k raw primers | `primers_filtered.{tsv,fasta}` | Tm+structure filter |
| primer3-py (multi) | Generate & Pre-filter | `test5.fasta` | `primer3py_output.txt` | Per-genome design |
| MELTING | Filter & Score | 50 primer seqs | `melting_results.txt` | ΔH/ΔS/Tm (NN model) |
| oligo-melting | Filter & Score | 1,061 filtered primers | `oligomelting_results.txt` | Tm at 1M + 50mM Na⁺ |
| blastn | Filter & Score | filtered primers vs 5 genomes | `blast_results.txt` | Off-target hits |
| MFEprimer | Filter & Score | filtered primers | `mfeprimer_{spec,dimer,hairpin}.txt` | Specificity + structure |
| seqkit | Filter & Score | filtered FASTA | `seqkit_primer_stats.txt` | FASTA stats |
| bowtie2 | Filter & Score | filtered primers vs 5 genomes | `bowtie2_primers.sam` | Alignment specificity |
| isPcr | Filter & Score | top 100 pairs | `ispcr_output.fasta` | Exact amplicons |
| ipcress | Filter & Score | top 100 pairs | `ipcress_output.txt` | Mismatch-tolerant PCR |
| primersearch | Filter & Score | top 100 pairs | `primersearch_output.txt` | EMBOSS PCR (10% mm) |
| tntblast | Filter & Score | 1 pair | `tntblast_output.txt` | Thermodynamic PCR sim |
| olivar | Optimize / Panel Selection | 3kb FASTA | `olivar_out/sars2_tiling.csv` + BED | SADDLE tiling |
| NGS-PrimerPlex | Optimize / Panel Selection | BED targets + genome | `*_combination_1.fa` | Multiplex panel |
| PrimalScheme3 | Specialized Workflows | MSA | `ps3_out/primer.bed` + `amplicon.bed` | Tiling scheme |
| varvamp | Specialized Workflows | MSA | `varvamp_out/primers.tsv` + BED | Variation-aware tiling |
| PUPpy | Specialized Workflows | CDS target + nontarget dirs | `ResultDB.tsv` → `UniquePrimerTable.tsv` | Taxon-specific |
| DegePrime | Specialized Workflows | MSA | `degeprime_output.tsv` | Degenerate primers |
| PrimerServer2 | Specialized Workflows | bracket-format FASTA | TSV primer pairs per region | Target-focused design |
| primerdiffer | Specialized Workflows | genome pairs | `primerdiffer_all/` | Discriminatory primers |
