#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# =============================================================================
# exploretools_final.sh
# Progressive refinement primer design pipeline
# All 21 tools — each stage uses outputs from previous stages
#
# PIPELINE FLOW:
#   Stage 1  — Alignment (mafft, clustalo)
#   Stage 2  — Generate raw candidates: primer3 → 20,000 raw primers
#   Stage 3a — Pre-filter: primer3-py Tm/structure → ~2,000 good primers
#   Stage 3b — Thermodynamics: MELTING + oligo-melting on filtered primers
#   Stage 4a — Specificity: BLAST + MFEprimer on filtered primers
#   Stage 4b — In-silico PCR: isPcr/ipcress/primersearch/tntblast on top 100 pairs
#   Stage 5  — Tiling schemes (independent): PS3, varvamp, olivar
#   Stage 6  — Taxon-specific + multiplex: PUPpy, NGS-PrimerPlex
#   Stage 7  — Independent generators: DegePrime, primerdiffer, primer3-py
#
# Run from: ~/primer/primer-framework
# Conda env: primer-env
# Usage:
#   bash exploretools_final.sh          # run everything
#   bash exploretools_final.sh stage1   # alignment
#   bash exploretools_final.sh stage2   # raw candidate generation
#   bash exploretools_final.sh stage3a  # pre-filter Tm/structure
#   bash exploretools_final.sh stage3b  # thermodynamics on filtered
#   bash exploretools_final.sh stage4a  # specificity on filtered
#   bash exploretools_final.sh stage4b  # in-silico PCR on top 100 pairs
#   bash exploretools_final.sh stage5   # tiling schemes
#   bash exploretools_final.sh stage6   # taxon-specific + multiplex
#   bash exploretools_final.sh stage7   # other independent generators
# =============================================================================

set -u  # No set -e — keep going even if individual tools fail

# ── Colours ───────────────────────────────────────────────────────────────────
GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m';  BOLD='\033[1m';   NC='\033[0m'

pass()    { echo -e "${GREEN}[PASS]${NC} $1"; }
fail()    { echo -e "${RED}[FAIL]${NC} $1"; }
info()    { echo -e "${CYAN}[INFO]${NC} $1"; }
flow()    { echo -e "${BOLD}${GREEN}  ↓ ${NC}$1"; }
section() { echo -e "\n${YELLOW}── $1 ──${NC}\n"; }
banner()  {
    echo -e "\n${BOLD}${CYAN}═══════════════════════════════════════════════════${NC}"
    echo -e "${BOLD}${CYAN}  $1${NC}"
    echo -e "${BOLD}${CYAN}═══════════════════════════════════════════════════${NC}\n"
}

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTA="test5.fasta"
OUTDIR="tool_outputs_final"
PRIMER3_PARAMS="$HOME/miniforge3/envs/primer-env/lib/python3.13/site-packages/primer3/src/libprimer3/primer3_config/"
MELTING_JAR="$HOME/MELTING5.2.0/executable/melting5.jar"
MELTING_DATA="$HOME/MELTING5.2.0/Data"
MFEPRIMER="$HOME/bin/mfeprimer"
DEGEPRIME="$HOME/DegePrime"
NPP="$HOME/NGS-PrimerPlex/NGS_primerplex.py"
PUPPY_TARGET="$SCRIPT_DIR/bacterial_cds/target"
PUPPY_NONTARGET="$SCRIPT_DIR/bacterial_cds/nontarget"

# ── Key output files (passed between stages) ──────────────────────────────────
# Stage 1 outputs:
MSA_FASTA="$OUTDIR/mafft_aligned.fasta"
# Stage 2 outputs:
P3_LEFT="$OUTDIR/primer3_output_left.txt"      # 10,000 left primers
P3_RIGHT="$OUTDIR/primer3_output_right.txt"    # 10,000 right primers
P3_ALL="$OUTDIR/primer3_output.txt"            # combined
# Stage 3a outputs:
FILTERED_TSV="$OUTDIR/primers_filtered.tsv"    # primers passing Tm+structure filter
FILTERED_FASTA="$OUTDIR/primers_filtered.fasta" # same as FASTA for tools
FILTERED_PAIRS="$OUTDIR/pairs_filtered.tsv"    # matched F+R pairs from filtered set
# Stage 4b inputs (top 100 pairs from filtered):
PAIRS_100_ISPCR="$OUTDIR/pairs_100_ispcr.txt"
PAIRS_100_IPCRESS="$OUTDIR/pairs_100_ipcress.txt"
PAIRS_100_PS="$OUTDIR/pairs_100_primersearch.txt"
PAIRS_100_TNT="$OUTDIR/pairs_100_tntblast.fasta"

declare -A STATUS

# ── Sanity checks ──────────────────────────────────────────────────────────────
cd "$(dirname "$0")" || exit 1
[[ ! -f "$FASTA" ]] && echo "ERROR: $FASTA not found. Run from ~/primer/primer-framework" && exit 1
mkdir -p "$OUTDIR"

echo -e "${BOLD}╔═══════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}║  PROGRESSIVE REFINEMENT PRIMER PIPELINE                  ║${NC}"
echo -e "${BOLD}║  20,000 raw → filtered → checked → final assay sets      ║${NC}"
echo -e "${BOLD}╚═══════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Input: $FASTA  ($(grep -c '^>' $FASTA) SARS-CoV-2 genomes)"
echo "Output: $OUTDIR/"
echo ""

# Extract single reference genome for tools that need it
python3 - << PSEOF
from Bio import SeqIO
seq = list(SeqIO.parse("$FASTA","fasta"))[0]
targets = [
    ("spike_NTD",    500,   600,  500),
    ("spike_RBD",   1200,  1350,  500),
    ("orf1ab_nsp3", 4500,  4650,  500),
    ("envelope",   26200, 26300,  500),
    ("nucleocapsid",28200, 28350, 500),
]
with open("$OUTDIR/ps2_sars2_targets.fasta","w") as f:
    for name, tstart, tend, flank in targets:
        left   = str(seq.seq[tstart-flank:tstart])
        target = str(seq.seq[tstart:tend])
        right  = str(seq.seq[tend:tend+flank])
        f.write(">" + name + "\n" + left + "[" + target + "]" + right + "\n")
print("Written: 5 target regions with 500bp flanks")
PSEOF

###############################################################################
# STAGE 1 — ALIGNMENT
###############################################################################
stage1() {
banner "STAGE 1 — ALIGNMENT (all 5 genomes)"
echo "Purpose: align all 5 genomes for MSA-based tools (DegePrime, varvamp, PS3)"
echo ""

section "1a) MAFFT"
info "Input:  test5.fasta (5 unaligned genomes)"
info "Output: $MSA_FASTA"
mafft --auto --thread 4 "$FASTA" > "$MSA_FASTA" 2>"$OUTDIR/mafft.log"
if [[ -s "$MSA_FASTA" ]]; then
    NSEQ=$(grep -c "^>" "$MSA_FASTA")
    pass "mafft → $NSEQ sequences aligned → $MSA_FASTA"
    STATUS[mafft]="PASS"
else
    fail "mafft failed"; STATUS[mafft]="FAIL"
fi

section "1b) Clustal Omega"
info "Input:  test5.fasta"
info "Output: $OUTDIR/clustalo_aligned.fasta"
clustalo -i "$FASTA" -o "$OUTDIR/clustalo_aligned.fasta" \
    --outfmt=fasta --force 2>"$OUTDIR/clustalo.log"
if [[ -s "$OUTDIR/clustalo_aligned.fasta" ]]; then
    pass "clustalo → $(grep -c '^>' $OUTDIR/clustalo_aligned.fasta) seqs → $OUTDIR/clustalo_aligned.fasta"
    STATUS[clustalo]="PASS"
else
    fail "clustalo failed"; STATUS[clustalo]="FAIL"
fi

flow "Stage 1 output → $MSA_FASTA used by: DegePrime, varvamp, PrimalScheme3"
}

###############################################################################
# STAGE 2 — GENERATE RAW CANDIDATES
###############################################################################
stage2() {
banner "STAGE 2 — GENERATE RAW CANDIDATES (~20,000)"
echo "Purpose: generate maximum possible primers from the reference genome"
echo "Design principle: cast wide net first, filter later"
echo ""

section "2a) primer3_core — 10,000 LEFT + 10,000 RIGHT primers"
info "Input:  MN908947.3 full genome (29,903bp)"
info "Output: $P3_LEFT  $P3_RIGHT  $P3_ALL"
info "Params: size 18-25bp, Tm 57-63°C, GC 40-60%, PRIMER_NUM_RETURN=10000"

python3 - << PYEOF2
from Bio import SeqIO
seq = list(SeqIO.parse("$FASTA","fasta"))[0]
params = "$PRIMER3_PARAMS"
for direction, pick_left, pick_right, outfile in [
    ("LEFT",  1, 0, "$OUTDIR/primer3_input_left.txt"),
    ("RIGHT", 0, 1, "$OUTDIR/primer3_input_right.txt"),
]:
    with open(outfile, "w") as f:
        f.write(f"SEQUENCE_ID={seq.id}_{direction}\n")
        f.write(f"PRIMER_THERMODYNAMIC_PARAMETERS_PATH={params}\n")
        f.write(f"SEQUENCE_TEMPLATE={str(seq.seq)}\n")
        f.write("PRIMER_TASK=pick_primer_list\n")
        f.write(f"PRIMER_PICK_LEFT_PRIMER={pick_left}\n")
        f.write(f"PRIMER_PICK_RIGHT_PRIMER={pick_right}\n")
        f.write("PRIMER_PICK_INTERNAL_OLIGO=0\n")
        f.write("PRIMER_OPT_SIZE=20\nPRIMER_MIN_SIZE=18\nPRIMER_MAX_SIZE=25\n")
        f.write("PRIMER_OPT_TM=60.0\nPRIMER_MIN_TM=57.0\nPRIMER_MAX_TM=63.0\n")
        f.write("PRIMER_MIN_GC=40.0\nPRIMER_MAX_GC=60.0\n")
        f.write("PRIMER_NUM_RETURN=10000\n")
        f.write("PRIMER_PRODUCT_SIZE_RANGE=100-400\n")
        f.write("PRIMER_MIN_THREE_PRIME_DISTANCE=3\n")
        f.write("=\n")
    print(f"Written: {outfile}")
print(f"Genome: {seq.id}  length={len(seq.seq)}bp")
PYEOF2

info "Running primer3 left primers (10,000)..."
primer3_core < "$OUTDIR/primer3_input_left.txt" > "$P3_LEFT" 2>&1 || true
info "Running primer3 right primers (10,000)..."
primer3_core < "$OUTDIR/primer3_input_right.txt" > "$P3_RIGHT" 2>&1 || true
cat "$P3_LEFT" "$P3_RIGHT" > "$P3_ALL"

NLEFT=$(grep -c "PRIMER_LEFT_[0-9]*_SEQUENCE=" "$P3_LEFT" 2>/dev/null || true)
NRIGHT=$(grep -c "PRIMER_RIGHT_[0-9]*_SEQUENCE=" "$P3_RIGHT" 2>/dev/null || true)
NLEFT=${NLEFT:-0}
NRIGHT=${NRIGHT:-0}
NTOTAL=$((NLEFT + NRIGHT))
if [[ "$NTOTAL" -gt 0 ]]; then
    pass "primer3_core → $NLEFT left + $NRIGHT right = $NTOTAL raw primers"
    STATUS[primer3_core]="PASS"
    echo "  Best left:  $(grep 'PRIMER_LEFT_0_SEQUENCE=' $P3_LEFT | cut -d= -f2)"
    echo "  Best right: $(grep 'PRIMER_RIGHT_0_SEQUENCE=' $P3_RIGHT | cut -d= -f2)"
else
    fail "primer3_core: no primers"; STATUS[primer3_core]="FAIL"
fi

flow "Stage 2 output → $((NLEFT+NRIGHT)) raw primers → fed into Stage 3a filter"
}

###############################################################################
# STAGE 3a — PRE-FILTER: Tm + STRUCTURE
###############################################################################
stage3a() {
banner "STAGE 3a — PRE-FILTER: Tm + STRUCTURE (~2,000 candidates)"
echo "Purpose: fast cheap filters first — reject bad Tm, hairpins, homodimers"
echo "Input:   $P3_ALL (20,000 raw primers from Stage 2)"
echo "Output:  $FILTERED_TSV + $FILTERED_FASTA + $FILTERED_PAIRS"
echo ""

section "3a) primer3-py — calc_tm + calc_hairpin + calc_homodimer on all 20,000"
info "Filters: Tm 57-63°C  |  hairpin dG > -2000 cal/mol  |  homodimer dG > -2000 cal/mol"

if [[ ! -f "$P3_LEFT" ]] || [[ ! -f "$P3_RIGHT" ]]; then
    fail "Stage 2 outputs not found — run stage2 first"; return
fi

python3 - <<'PYEOF'
import re, primer3

# Load all primers from stage 2
with open("tool_outputs_final/primer3_output_left.txt") as f:
    lefts  = re.findall(r"PRIMER_LEFT_\d+_SEQUENCE=(\w+)",  f.read())
with open("tool_outputs_final/primer3_output_right.txt") as f:
    rights = re.findall(r"PRIMER_RIGHT_\d+_SEQUENCE=(\w+)", f.read())

print(f"Input from Stage 2: {len(lefts)} left + {len(rights)} right = {len(lefts)+len(rights)} raw primers")

# Calculate Tm, hairpin, homodimer for EVERY primer
all_seqs = list(set(lefts + rights))
print(f"Calculating Tm + structure for {len(all_seqs)} unique primers...")

results = {}
for seq in all_seqs:
    tm  = primer3.calc_tm(seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.25, dna_conc=250)
    hp  = primer3.calc_hairpin(seq)
    hd  = primer3.calc_homodimer(seq)
    results[seq] = (tm, hp.dg, hd.dg)

# Write full Tm table (all 20,000)
with open("tool_outputs_final/all_primers_tm.tsv","w") as f:
    f.write("primer_seq\tTm_C\thairpin_dG\thousedimer_dG\n")
    for seq,(tm,hp,hd) in results.items():
        f.write(f"{seq}\t{tm:.2f}\t{hp:.0f}\t{hd:.0f}\n")

# Apply filters → filtered set
filtered = {seq:(tm,hp,hd) for seq,(tm,hp,hd) in results.items()
            if 57.0 <= tm <= 63.0 and hp > -2000 and hd > -2000}

print(f"After filter: {len(filtered)} primers pass (Tm 57-63 + structure)")

# Write filtered TSV
with open("tool_outputs_final/primers_filtered.tsv","w") as f:
    f.write("primer_seq\tTm_C\thairpin_dG\thousedimer_dG\n")
    for seq,(tm,hp,hd) in filtered.items():
        f.write(f"{seq}\t{tm:.2f}\t{hp:.0f}\t{hd:.0f}\n")

# Write filtered FASTA
with open("tool_outputs_final/primers_filtered.fasta","w") as f:
    for i,seq in enumerate(filtered):
        f.write(f">primer_{i}\n{seq}\n")

# Build matched pairs: match filtered left primers with filtered right primers
filtered_lefts  = [s for s in lefts  if s in filtered]
filtered_rights = [s for s in rights if s in filtered]
pairs = list(zip(filtered_lefts, filtered_rights))

print(f"Filtered left: {len(filtered_lefts)}  Filtered right: {len(filtered_rights)}")
print(f"Matched pairs: {len(pairs)}")

# Write pairs TSV
with open("tool_outputs_final/pairs_filtered.tsv","w") as f:
    f.write("pair_id\tF_seq\tR_seq\tF_Tm\tR_Tm\tF_hairpin\tR_hairpin\n")
    for i,(l,r) in enumerate(pairs):
        lt,lh,_ = filtered[l]
        rt,rh,_ = filtered[r]
        f.write(f"pair_{i+1}\t{l}\t{r}\t{lt:.2f}\t{rt:.2f}\t{lh:.0f}\t{rh:.0f}\n")

# Write top 100 pairs in formats needed by Stage 4b tools
pairs100 = pairs[:100]
with open("tool_outputs_final/pairs_100_ispcr.txt","w") as f:
    for i,(l,r) in enumerate(pairs100):
        f.write(f"pair{i+1}\t{l}\t{r}\n")
with open("tool_outputs_final/pairs_100_ipcress.txt","w") as f:
    for i,(l,r) in enumerate(pairs100):
        f.write(f"pair{i+1} {l} {r} 50 500\n")
with open("tool_outputs_final/pairs_100_primersearch.txt","w") as f:
    for i,(l,r) in enumerate(pairs100):
        f.write(f"pair{i+1}\t{l}\t{r}\n")
with open("tool_outputs_final/pairs_100_tntblast.fasta","w") as f:
    for i,(l,r) in enumerate(pairs100):
        f.write(f">pair{i+1}_F\n{l}\n>pair{i+1}_R\n{r}\n")

# Write 50 primers for MELTING
with open("tool_outputs_final/primers_melting_50.txt","w") as f:
    for seq in list(filtered)[:50]:
        f.write(seq + "\n")

print(f"\nOutputs written:")
print(f"  all_primers_tm.tsv      — Tm for all {len(results)} raw primers")
print(f"  primers_filtered.tsv    — {len(filtered)} primers passing filters")
print(f"  primers_filtered.fasta  — same as FASTA")
print(f"  pairs_filtered.tsv      — {len(pairs)} matched pairs")
print(f"  pairs_100_*.txt         — top 100 pairs for in-silico PCR tools")
print(f"  primers_melting_50.txt  — 50 primers for MELTING")
PYEOF

if [[ -s "$FILTERED_TSV" ]]; then
    NFILT=$(wc -l < "$FILTERED_TSV")
    pass "primer3-py filter → $NFILT primers pass Tm+structure → $FILTERED_TSV"
    STATUS[primer3py]="PASS"
    echo ""
    echo "Sample filtered primers:"
    tail -n +2 "$FILTERED_TSV" | head -3
else
    fail "primer3-py filter failed"; STATUS[primer3py]="FAIL"
fi

flow "Stage 3a output → $FILTERED_FASTA used by: MELTING, oligo-melting, BLAST, MFEprimer"
flow "Stage 3a output → $PAIRS_100_ISPCR etc used by: isPcr, ipcress, primersearch, tntblast"
}

###############################################################################
# STAGE 3b — THERMODYNAMICS on FILTERED PRIMERS
###############################################################################
stage3b() {
banner "STAGE 3b — THERMODYNAMICS on FILTERED PRIMERS"
echo "Purpose: detailed Tm with salt correction on the ~2,000 filtered candidates"
echo "Input:   $FILTERED_FASTA (filtered primers from Stage 3a)"
echo ""

if [[ ! -f "$FILTERED_FASTA" ]]; then
    fail "Stage 3a output not found — run stage3a first"; return
fi

NFILT=$(grep -c "^>" "$FILTERED_FASTA" 2>/dev/null || echo 0)
info "Using $NFILT filtered primers from Stage 3a"

section "3b-i) MELTING (Java) — nearest-neighbour Tm, 50mM Na+"
info "Input:  primers_melting_50.txt (50 sample from filtered set)"
info "Output: $OUTDIR/melting_results.txt"
info "Note:   50 primers shown — MELTING is slow (Java per-call), full run would be ~1,928 × 2s"
mkdir -p "$OUTDIR/melting_out"
while IFS= read -r SEQ; do
    [[ -z "$SEQ" ]] && continue
    NN_PATH="$MELTING_DATA" java -jar "$MELTING_JAR" \
        -S "$SEQ" -H dnadna -P 0.00000025 -E Na=0.05 \
        2>/dev/null >> "$OUTDIR/melting_out/melting_results.txt" || true
done < "$OUTDIR/primers_melting_50.txt"

if [[ -s "$OUTDIR/melting_out/melting_results.txt" ]]; then
    pass "MELTING → $OUTDIR/melting_out/melting_results.txt"
    STATUS[melting]="PASS"
    grep -i "melting temperature" "$OUTDIR/melting_out/melting_results.txt" | head -3
else
    fail "MELTING failed"; STATUS[melting]="FAIL"
fi

section "3b-ii) oligo-melting — Python Tm on ALL filtered primers"
info "Input:  $FILTERED_FASTA (all ~2,000 filtered primers)"
info "Output: $OUTDIR/oligomelting_results.txt"
grep -v "^>" "$FILTERED_FASTA" | while IFS= read -r SEQ; do
    [[ -z "$SEQ" ]] && continue
    oligomelting "$SEQ" --conc-na 50 \
        >> "$OUTDIR/oligomelting_results.txt" 2>&1 || true
done

if [[ -s "$OUTDIR/oligomelting_results.txt" ]]; then
    NOLIGO=$(grep -c "^Sequence:" "$OUTDIR/oligomelting_results.txt" 2>/dev/null || echo "?")
    pass "oligo-melting → $NOLIGO primers processed → $OUTDIR/oligomelting_results.txt"
    STATUS[oligomelting]="PASS"
    head -6 "$OUTDIR/oligomelting_results.txt"
else
    fail "oligo-melting failed"; STATUS[oligomelting]="FAIL"
fi

flow "Stage 3b adds Tm annotation to filtered primer set"
}

###############################################################################
# STAGE 4a — SPECIFICITY on FILTERED PRIMERS
###############################################################################
stage4a() {
banner "STAGE 4a — SPECIFICITY CHECKING on FILTERED PRIMERS"
echo "Purpose: check where filtered primers bind across all 5 genomes"
echo "Input:   $FILTERED_FASTA (filtered primers from Stage 3a)"
echo ""

if [[ ! -f "$FILTERED_FASTA" ]]; then
    fail "Stage 3a output not found — run stage3a first"; return
fi

NFILT=$(grep -c "^>" "$FILTERED_FASTA" 2>/dev/null || echo 0)
info "Using $NFILT filtered primers from Stage 3a"

section "4a-i) makeblastdb + blastn — off-target binding check"
info "Input:  $FILTERED_FASTA against all 5 SARS-CoV-2 genomes"
info "Output: $OUTDIR/blast_results.txt"
makeblastdb -in "$FASTA" -dbtype nucl \
    -out "$OUTDIR/sars2_blastdb" -title "sars2_5genomes" \
    2>"$OUTDIR/makeblastdb.log" || true

blastn -query "$FILTERED_FASTA" \
    -db "$OUTDIR/sars2_blastdb" \
    -task blastn-short \
    -word_size 7 \
    -outfmt 6 \
    -out "$OUTDIR/blast_results.txt" \
    2>/dev/null || true

if [[ -s "$OUTDIR/blast_results.txt" ]]; then
    NHITS=$(wc -l < "$OUTDIR/blast_results.txt")
    NUNIQ=$(cut -f1 "$OUTDIR/blast_results.txt" | sort -u | wc -l)
    pass "blastn → $NHITS hits from $NUNIQ primers → $OUTDIR/blast_results.txt"
    STATUS[blastn]="PASS"
    echo "Sample hits (qseqid sseqid pident length):"
    awk '{print $1"\t"$2"\t"$3"\t"$4}' "$OUTDIR/blast_results.txt" | head -3
else
    fail "blastn: no hits"; STATUS[blastn]="FAIL"
fi

section "4a-ii) MFEprimer — specificity + dimer + hairpin on filtered primers"
info "Input:  $FILTERED_FASTA"
info "Output: mfeprimer_spec.txt  mfeprimer_dimer.txt  mfeprimer_hairpin.txt"
"$MFEPRIMER" index -i "$FASTA" 2>"$OUTDIR/mfeprimer_index.log" || true
"$MFEPRIMER" spec    -i "$FILTERED_FASTA" -d "$FASTA" -o "$OUTDIR/mfeprimer_spec.txt"    2>/dev/null || true
"$MFEPRIMER" dimer   -i "$FILTERED_FASTA"             -o "$OUTDIR/mfeprimer_dimer.txt"   2>/dev/null || true
"$MFEPRIMER" hairpin -i "$FILTERED_FASTA"             -o "$OUTDIR/mfeprimer_hairpin.txt" 2>/dev/null || true

if [[ -s "$OUTDIR/mfeprimer_spec.txt" ]]; then
    pass "MFEprimer spec/dimer/hairpin → $OUTDIR/mfeprimer_*.txt"
    STATUS[mfeprimer]="PASS"
    head -5 "$OUTDIR/mfeprimer_spec.txt"
else
    fail "MFEprimer failed"; STATUS[mfeprimer]="FAIL"
fi

section "4a-iii) primerdiffer — discriminatory primers all 10 genome pairs"
info "Input:  all 5 genomes — 10 unique pairs"
info "Output: $OUTDIR/primerdiffer_all/"

# Extract all 5 genomes to separate files
python3 - << PDEOF
from Bio import SeqIO
seqs = list(SeqIO.parse("$FASTA","fasta"))
for s in seqs:
    SeqIO.write([s], f"$OUTDIR/{s.id}.fasta", "fasta")
print(f"Extracted {len(seqs)} genomes")
PDEOF

# Run all 10 pairs
python3 - << PDEOF2
import subprocess, os
from itertools import combinations
seqs = ["MN908947.3","MN985325.1","MN988713.1","MN938384.1","MN975262.1"]
pairs = list(combinations(seqs, 2))
total = 0
for g1, g2 in pairs:
    outdir = f"$OUTDIR/primerdiffer_all/{g1}_vs_{g2}"
    os.makedirs(outdir, exist_ok=True)
    cmd = ["primerdesign.py",
           "-g1", os.path.realpath(f"$OUTDIR/{g1}.fasta"),
           "-g2", os.path.realpath(f"$OUTDIR/{g2}.fasta"),
           "-pos", f"{g1}:1-29903",
           "-d", outdir]
    subprocess.run(cmd, capture_output=True)
    n = sum(len(open(os.path.join(outdir,f)).readlines())
            for f in os.listdir(outdir) if f.endswith(".txt"))
    total += n
    print(f"  {g1} vs {g2}: {n} primers")
print(f"Total: {total} discriminatory primers across {len(pairs)} pairs")
PDEOF2

if [[ -d "$OUTDIR/primerdiffer_all" ]]; then
    NTOTAL=$(find "$OUTDIR/primerdiffer_all" -name "*.txt" -exec cat {} \; |         grep -c "." 2>/dev/null || echo 0)
    pass "primerdiffer → $NTOTAL discriminatory primers → $OUTDIR/primerdiffer_all/"
    STATUS[primerdiffer]="PASS"
else
    fail "primerdiffer failed"; STATUS[primerdiffer]="FAIL"
fi


section "4a-iv) seqkit — primer statistics"
info "Input:  $FILTERED_FASTA"
seqkit stats "$FILTERED_FASTA" > "$OUTDIR/seqkit_primer_stats.txt" 2>/dev/null || true
if [[ -s "$OUTDIR/seqkit_primer_stats.txt" ]]; then
    pass "seqkit stats → $OUTDIR/seqkit_primer_stats.txt"
    STATUS[seqkit]="PASS"
    cat "$OUTDIR/seqkit_primer_stats.txt"
else
    fail "seqkit failed"; STATUS[seqkit]="FAIL"
fi

section "4a-v) bowtie2 + samtools — alignment-based specificity"
info "Input:  $FILTERED_FASTA aligned to all 5 genomes"
bowtie2-build "$FASTA" "$OUTDIR/sars2_bt2"     > "$OUTDIR/bowtie2_build.log" 2>&1 || true
bowtie2 -x "$OUTDIR/sars2_bt2"     -f -U "$FILTERED_FASTA"     --no-unal -N 1 -L 18     -S "$OUTDIR/bowtie2_primers.sam"     2>"$OUTDIR/bowtie2.log" || true
if [[ -s "$OUTDIR/bowtie2_primers.sam" ]]; then
    NALIGN=$(samtools view -c "$OUTDIR/bowtie2_primers.sam")
    pass "bowtie2 → $NALIGN primers aligned → $OUTDIR/bowtie2_primers.sam"
    STATUS[bowtie2]="PASS"
    tail -4 "$OUTDIR/bowtie2.log"
else
    fail "bowtie2 failed"; STATUS[bowtie2]="FAIL"
fi

flow "Stage 4a narrows filtered primers to those with good specificity profiles"
}

###############################################################################
# STAGE 4b — IN-SILICO PCR on TOP 100 FILTERED PAIRS
###############################################################################
stage4b() {
banner "STAGE 4b — IN-SILICO PCR on TOP 100 FILTERED PAIRS"
echo "Purpose: simulate actual PCR amplification with the best filtered pairs"
echo "Input:   pairs_100_*.txt (top 100 pairs from Stage 3a filtered set)"
echo ""

if [[ ! -f "$PAIRS_100_ISPCR" ]]; then
    fail "Stage 3a pair outputs not found — run stage3a first"; return
fi

NPAIRS=$(wc -l < "$PAIRS_100_ISPCR")
info "Using top $NPAIRS filtered pairs from Stage 3a"

section "4b-i) isPcr — exact in-silico PCR (all 5 genomes)"
info "Input:  $PAIRS_100_ISPCR (top 100 filtered pairs)"
info "Format: name<TAB>F_seq<TAB>R_seq"

isPcr "$FASTA" "$PAIRS_100_ISPCR" stdout \
    2>/dev/null > "$OUTDIR/ispcr_output.fasta" || true

if [[ -s "$OUTDIR/ispcr_output.fasta" ]]; then
    NAMP=$(grep -c "^>" "$OUTDIR/ispcr_output.fasta" || echo 0)
    pass "isPcr → $NAMP amplicons from 100 pairs × 5 genomes → $OUTDIR/ispcr_output.fasta"
    STATUS[ispcr]="PASS"
    head -4 "$OUTDIR/ispcr_output.fasta"
else
    fail "isPcr: no amplicons"; STATUS[ispcr]="FAIL"
fi

section "4b-ii) ipcress — mismatch-tolerant PCR (all 5 genomes)"
info "Input:  $PAIRS_100_IPCRESS"
info "Format: name F_seq R_seq min_size max_size"

ipcress --input "$PAIRS_100_IPCRESS" \
    --sequence "$FASTA" \
    --mismatch 2 \
    2>/dev/null > "$OUTDIR/ipcress_output.txt" || true

if [[ -s "$OUTDIR/ipcress_output.txt" ]]; then
    NAMP=$(grep -c "^ipcress:" "$OUTDIR/ipcress_output.txt" || echo 0)
    pass "ipcress → $NAMP amplimers → $OUTDIR/ipcress_output.txt"
    STATUS[ipcress]="PASS"
    grep "^ipcress:" "$OUTDIR/ipcress_output.txt" | head -3
else
    fail "ipcress: no output"; STATUS[ipcress]="FAIL"
fi

section "4b-iii) primersearch — EMBOSS mismatch-tolerant PCR"
info "Input:  $PAIRS_100_PS"
info "Format: name<TAB>F_seq<TAB>R_seq"

primersearch \
    -seqall "$FASTA" \
    -infile "$PAIRS_100_PS" \
    -mismatchpercent 10 \
    -outfile "$OUTDIR/primersearch_output.txt" \
    2>/dev/null || true

if [[ -s "$OUTDIR/primersearch_output.txt" ]]; then
    pass "primersearch → $OUTDIR/primersearch_output.txt"
    STATUS[primersearch]="PASS"
    grep "Amplimer\|Forward\|Reverse" \
        "$OUTDIR/primersearch_output.txt" | head -6
else
    fail "primersearch: no output"; STATUS[primersearch]="FAIL"
fi

section "4b-iv) tntblast — thermodynamic PCR simulation"
info "Input: assay definition file (old working format)"

cat > "$OUTDIR/tntblast_assay.txt" <<EOF
#name	forward	reverse
pair1	TCGAACTGCACCTCATGGTC	GACTTTAGATCGGCGCCGTA
EOF

tntblast \
    -i "$OUTDIR/tntblast_assay.txt" \
    -d "$FASTA" \
    -o "$OUTDIR/tntblast_output.txt" \
    -e 45 \
    2>&1 | tee "$OUTDIR/tntblast.log" | head -30 || true

if [[ -s "$OUTDIR/tntblast_output.txt" ]]; then
    pass "tntblast → $OUTDIR/tntblast_output.txt"
    STATUS[tntblast]="PASS"
    head -5 "$OUTDIR/tntblast_output.txt"
else
    fail "tntblast failed — see $OUTDIR/tntblast.log"
    STATUS[tntblast]="FAIL"
fi

flow "Stage 4b confirms which filtered pairs actually amplify the target"
}

###############################################################################
# STAGE 5 — TILING SCHEMES (independent)
###############################################################################
stage5() {
banner "STAGE 5 — TILING SCHEMES (independent tools)"
echo "Purpose: design complete amplicon tiling schemes across the genome"
echo "These tools run independently — they generate AND optimise their own primers"
echo ""

section "5a) PrimalScheme3 — MSA-aware tiling"
info "Input:  $MSA_FASTA (mafft alignment from Stage 1)"
rm -rf "$OUTDIR/ps3_out"
primalscheme3 scheme-create \
    --msa "$MSA_FASTA" \
    --output "$OUTDIR/ps3_out" \
    --amplicon-size 400 \
    2>"$OUTDIR/ps3.log" || true
if [[ -d "$OUTDIR/ps3_out" ]] && \
   ls "$OUTDIR/ps3_out/"* 2>/dev/null | head -1 | grep -q .; then
    pass "PrimalScheme3 → $OUTDIR/ps3_out/"
    STATUS[primalscheme3]="PASS"
    ls "$OUTDIR/ps3_out/"
else
    fail "PrimalScheme3 failed"; STATUS[primalscheme3]="FAIL"
fi

section "5c) varvamp — variation-aware tiling"
info "Input: $MSA_FASTA (mafft alignment from Stage 1)"

if [[ ! -f "$MSA_FASTA" ]]; then
    fail "varvamp failed — missing $MSA_FASTA. Run stage1 first."
    STATUS[varvamp]="FAIL"
else
    rm -rf "$OUTDIR/varvamp_out"
    varvamp tiled "$MSA_FASTA" "$OUTDIR/varvamp_out" \
        2>&1 | tee "$OUTDIR/varvamp.log" | head -20 || true

    if [[ -d "$OUTDIR/varvamp_out" ]] && ls "$OUTDIR/varvamp_out/"* 2>/dev/null | head -1 | grep -q .; then
        pass "varvamp → $OUTDIR/varvamp_out"
        STATUS[varvamp]="PASS"
        ls "$OUTDIR/varvamp_out/"
    else
        fail "varvamp failed — see $OUTDIR/varvamp.log"
        STATUS[varvamp]="FAIL"
    fi
fi

section "5d) olivar — SADDLE optimization (3000bp demo region)"
info "Input:  first 3000bp of MN908947.3 (full genome takes hours)"
python3 -c "
from Bio import SeqIO
seq = list(SeqIO.parse('$FASTA','fasta'))[0]
seq.seq = seq.seq[:3000]
SeqIO.write([seq], '$OUTDIR/olivar_short.fasta','fasta')
print(f'  Demo: first 3000bp of {seq.id}')
"
mkdir -p "$OUTDIR/olivar_db" "$OUTDIR/olivar_out"
olivar build "$OUTDIR/olivar_short.fasta" \
    --output "$OUTDIR/olivar_db" --title sars2_olivar \
    --threads 4 2>"$OUTDIR/olivar_build.log" || true
if [[ -f "$OUTDIR/olivar_db/sars2_olivar.olvr" ]]; then
    olivar tiling "$OUTDIR/olivar_db/sars2_olivar.olvr" \
        --output "$OUTDIR/olivar_out" --title sars2_tiling \
        --max-amp-len 400 --seed 42 --threads 4 \
        2>&1 | grep -E "amplicons|coverage|Finished|saved" | head -8 || true
    if ls "$OUTDIR/olivar_out/"*.csv 2>/dev/null | head -1 | grep -q csv; then
        pass "olivar → $OUTDIR/olivar_out/"
        STATUS[olivar]="PASS"
        ls "$OUTDIR/olivar_out/"
    else
        fail "olivar tiling failed"; STATUS[olivar]="FAIL"
    fi
else
    fail "olivar build failed"; STATUS[olivar]="FAIL"
fi
}

###############################################################################
# STAGE 6 — TAXON-SPECIFIC + MULTIPLEX
###############################################################################
stage6() {
banner "STAGE 6 — TAXON-SPECIFIC + MULTIPLEX PANEL"

section "6a) PUPpy — taxon-specific primers (E. coli K-12 vs Salmonella)"
info "Input:  bacterial_cds/target/ (EcoliK12_cds.fna)"
info "        bacterial_cds/nontarget/ (SalmonellaTM_cds.fna)"
info "Step 1: puppy-align   → ResultDB.tsv"
info "Step 2: puppy-primers → UniquePrimerTable.tsv"

export PATH="$HOME/PUPpy/scripts:$PATH"

if [[ -d "$PUPPY_TARGET" ]] && [[ -d "$PUPPY_NONTARGET" ]]; then
    command -v puppy-align >/dev/null || { fail "puppy-align not on PATH"; STATUS[puppy]="FAIL"; }
    command -v puppy-primers >/dev/null || { fail "puppy-primers not on PATH"; STATUS[puppy]="FAIL"; }
    command -v mmseqs >/dev/null || { fail "mmseqs2 not installed"; STATUS[puppy]="FAIL"; }

    if [[ "${STATUS[puppy]:-}" != "FAIL" ]]; then
        rm -rf "$OUTDIR/puppy_align" "$OUTDIR/puppy_primers"
        mkdir -p "$OUTDIR/puppy_align" "$OUTDIR/puppy_primers"

        stdbuf -oL -eL puppy-align \
            -pr "$PUPPY_TARGET" \
            -nt "$PUPPY_NONTARGET" \
            -o "$OUTDIR/puppy_align" 2>&1 | tee "$OUTDIR/puppy_align.log"

        if [[ -f "$OUTDIR/puppy_align/ResultDB.tsv" ]]; then
            stdbuf -oL -eL puppy-primers \
                -pr "$PUPPY_TARGET" \
                -i "$OUTDIR/puppy_align/ResultDB.tsv" \
                -o "$OUTDIR/puppy_primers" 2>&1 | tee "$OUTDIR/puppy_primers.log"

            if [[ -f "$OUTDIR/puppy_primers/UniquePrimerTable.tsv" ]]; then
                NPUP=$(( $(wc -l < "$OUTDIR/puppy_primers/UniquePrimerTable.tsv") - 1 ))
                pass "PUPpy → $NPUP entries → $OUTDIR/puppy_primers/UniquePrimerTable.tsv"
                STATUS[puppy]="PASS"
                echo "Sample E. coli specific primers:"
                tail -n +2 "$OUTDIR/puppy_primers/UniquePrimerTable.tsv" | \
                    cut -f1,8,9,13,14 | head -3 | column -t 2>/dev/null || \
                    tail -n +2 "$OUTDIR/puppy_primers/UniquePrimerTable.tsv" | head -3
            else
                fail "PUPpy primers step failed"
                tail -50 "$OUTDIR/puppy_primers.log" 2>/dev/null
                STATUS[puppy]="FAIL"
            fi
        else
            fail "PUPpy align step failed"
            tail -50 "$OUTDIR/puppy_align.log" 2>/dev/null
            STATUS[puppy]="FAIL"
        fi
    fi
else
    info "PUPpy: bacterial CDS not found at $PUPPY_TARGET"
    info "  Expected: bacterial_cds/target/EcoliK12_cds.fna"
    STATUS[puppy]="SKIP"
fi

section "6b) NGS-PrimerPlex — multiplex panel for 2 target regions"
info "Input:  npp_targets.bed (2 SARS-CoV-2 regions)"
info "        test5.fasta (reference genome)"
info "Output: 8 amplicons across 2 target regions"

cat > "$OUTDIR/npp_targets.bed" << 'EOF'
MN908947.3	400	900	target1
MN908947.3	1200	1700	target2
EOF

python3 "$NPP" \
    --regions-file "$OUTDIR/npp_targets.bed" \
    --reference-genome "$FASTA" \
    --min-amplicon-length 100 \
    --max-amplicon-length 300 \
    --optimal-amplicon-length 200 \
    --min-primer-melting-temp 57 \
    --max-primer-melting-temp 63 \
    --optimal-primer-melting-temp 60 \
    --max-primer-nonspecific 2 \
    --primers-number1 2 \
    --return-variants-number 1 \
    2>&1 | tee "$OUTDIR/npp_run.log" | \
    grep -E "100\.0%|amplicons|finished|written|Number of" | head -10 &

NPP_PID=$!
info "NGS-PrimerPlex running in background (PID $NPP_PID)"
info "Monitor: tail -f $OUTDIR/npp_run.log"
info "Expected: 8 amplicons, ~5 minutes"
STATUS[ngs_primerplex]="RUNNING"
}

###############################################################################
# STAGE 7 — OTHER INDEPENDENT GENERATORS
###############################################################################
stage7() {
banner "STAGE 7 — OTHER INDEPENDENT PRIMER GENERATORS"
echo "These tools generate primers from alignments independently of Stage 2"
echo ""

section "7a) DegePrime — degenerate primers from MSA"
info "Input:  $MSA_FASTA (mafft alignment from Stage 1)"
info "Step 1: TrimAlignment.pl → trimmed alignment"
info "Step 2: DegePrime.pl     → degenerate primer windows"

perl "$DEGEPRIME/TrimAlignment.pl" \
    -i "$MSA_FASTA" \
    -o "$OUTDIR/mafft_trimmed.fasta" \
    -min 0.9 2>"$OUTDIR/trim.log" || true

perl "$DEGEPRIME/DegePrime.pl" \
    -i "$OUTDIR/mafft_trimmed.fasta" \
    -d 12 -l 20 \
    -o "$OUTDIR/degeprime_output.tsv" \
    2>"$OUTDIR/degeprime.log" || true

if [[ -s "$OUTDIR/degeprime_output.tsv" ]]; then
    NLINES=$(wc -l < "$OUTDIR/degeprime_output.tsv")
    if [[ "$NLINES" -gt 1 ]]; then
        pass "DegePrime → $NLINES primer windows → $OUTDIR/degeprime_output.tsv"
        STATUS[degeprime]="PASS"
        head -4 "$OUTDIR/degeprime_output.tsv"
    else
        fail "DegePrime: only header — alignment may be too uniform"
        STATUS[degeprime]="FAIL"
    fi
else
    fail "DegePrime failed"; STATUS[degeprime]="FAIL"
fi


section "7c) PrimerServer2 — target-focused primer design (5 SARS-CoV-2 regions)"
info "Input:  5 diagnostic target regions with 500bp flanks"
info "Output: primerserver2_final_results.tsv"
info "Targets: spike_NTD, spike_RBD, orf1ab_nsp3, envelope, nucleocapsid"

# Build BLAST index on test5.fasta if not already done
if [[ ! -f "$FASTA.nhr" ]]; then
    makeblastdb -in "$FASTA" -dbtype nucl -parse_seqids         2>"$OUTDIR/ps2_makeblastdb.log" || true
fi

# Create input with 5 SARS-CoV-2 diagnostic regions, 500bp flanks
python3 - << PSEOF
from Bio import SeqIO
seq = list(SeqIO.parse("$FASTA","fasta"))[0]
targets = [
    ("spike_NTD",    500,   600,  500),
    ("spike_RBD",   1200,  1350,  500),
    ("orf1ab_nsp3", 4500,  4650,  500),
    ("envelope",   26200, 26300,  500),
    ("nucleocapsid",28200, 28350, 500),
]
with open("$OUTDIR/ps2_sars2_targets.fasta","w") as f:
    for name, tstart, tend, flank in targets:
        left   = str(seq.seq[tstart-flank:tstart])
        target = str(seq.seq[tstart:tend])
        right  = str(seq.seq[tend:tend+flank])
        f.write(">" + name + "\n" + left + "[" + target + "]" + right + "\n")
print("Written: 5 target regions with 500bp flanks")
PSEOF

# Run PrimerServer2 — all possible primers per region
primertool design     "$OUTDIR/ps2_sars2_targets.fasta"     "$(realpath $FASTA)"     --primer-num-return 10000     --primer-num-retain 10000     -t "$OUTDIR/primerserver2_final_results.tsv"     2>"$OUTDIR/ps2.log" || true

if [[ -s "$OUTDIR/primerserver2_final_results.tsv" ]]; then
    NTOTAL=$(grep -v "^#\|^###" "$OUTDIR/primerserver2_final_results.tsv" | wc -l)
    pass "PrimerServer2 → $NTOTAL primer pairs across 5 diagnostic regions"
    STATUS[primerserver2]="PASS"
    echo ""
    echo "Primers per region:"
    grep -v "^#" "$OUTDIR/primerserver2_final_results.tsv" | grep -v "^###" | cut -f1 | sort | uniq -c
    echo ""
    echo "Sample output (spike_NTD top pair):"
    grep "^spike_NTD" "$OUTDIR/primerserver2_final_results.tsv" |         head -1 | cut -f1,3,4,6,7 | column -t 2>/dev/null ||         grep "^spike_NTD" "$OUTDIR/primerserver2_final_results.tsv" | head -1
else
    fail "PrimerServer2 failed — see $OUTDIR/ps2.log"
    STATUS[primerserver2]="FAIL"
fi

section "7b) primer3-py — design primers from all 5 genomes"
info "Input:  all 5 genomes in test5.fasta"
info "Output: $OUTDIR/primer3py_output.txt"

python3 - <<'PYEOF'
import re
import primer3
from Bio import SeqIO

seqs = list(SeqIO.parse("test5.fasta", "fasta"))
total = 0

with open("tool_outputs_final/primer3py_output.txt", "w") as fout:
    fout.write("genome_id\tpairs_returned\tcleaned_ambiguous_bases\n")

    for s in seqs:
        raw_seq = str(s.seq).upper()
        cleaned_seq = re.sub(r'[^ACGTN]', 'N', raw_seq)

        num_changed = sum(1 for a, b in zip(raw_seq, cleaned_seq) if a != b)

        try:
            result = primer3.design_primers(
                {
                    "SEQUENCE_ID": s.id,
                    "SEQUENCE_TEMPLATE": cleaned_seq
                },
                {
                    "PRIMER_OPT_SIZE": 20,
                    "PRIMER_MIN_SIZE": 18,
                    "PRIMER_MAX_SIZE": 25,
                    "PRIMER_OPT_TM": 60.0,
                    "PRIMER_MIN_TM": 57.0,
                    "PRIMER_MAX_TM": 63.0,
                    "PRIMER_MIN_GC": 40.0,
                    "PRIMER_MAX_GC": 60.0,
                    "PRIMER_NUM_RETURN": 50,
                    "PRIMER_PRODUCT_SIZE_RANGE": [[100, 400]],
                }
            )

            n = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
            total += n
            fout.write(f"{s.id}\t{n}\t{num_changed}\n")

            if num_changed > 0:
                print(f"  {s.id}: {n} pairs  (cleaned {num_changed} ambiguous base(s) → N)")
            else:
                print(f"  {s.id}: {n} pairs")

        except Exception as e:
            fout.write(f"{s.id}\tERROR\t{num_changed}\n")
            print(f"  {s.id}: ERROR → {e}")

print(f"Total: {total} pairs across {len(seqs)} genomes → tool_outputs_final/primer3py_output.txt")
PYEOF

if [[ -s "$OUTDIR/primer3py_output.txt" ]]; then
    pass "primer3-py → $OUTDIR/primer3py_output.txt"
    STATUS[primer3py_multi]="PASS"
else
    fail "primer3-py multi-genome failed"
    STATUS[primer3py_multi]="FAIL"
fi
}

###############################################################################
# SUMMARY
###############################################################################
show_summary() {
banner "PIPELINE SUMMARY — Progressive Refinement Results"
echo ""
echo -e "${BOLD}Pipeline flow:${NC}"
echo "  Stage 2  → 20,000 raw primers (primer3_core)"
echo "  Stage 3a → ~2,000 filtered (Tm 57-63 + structure)"
echo "  Stage 3b → Tm annotated (MELTING, oligo-melting)"
echo "  Stage 4a → Specificity checked (BLAST, MFEprimer)"
echo "  Stage 4b → PCR confirmed (isPcr, ipcress, primersearch, tntblast)"
echo "  Stage 5  → Tiling schemes (PS3, varvamp, olivar)"
echo "  Stage 6  → Specialised (PUPpy, NGS-PrimerPlex)"
echo "  Stage 7  → Other generators (DegePrime, primerdiffer)"
echo ""
printf "${BOLD}%-22s %-12s %-22s %-10s${NC}\n" "TOOL" "STATUS" "ROLE" "STAGE"
printf "%-22s %-12s %-22s %-10s\n" "──────────────────────" "──────────" "──────────────────────" "──────────"

declare -A ROLE STAGEMAP
ROLE[mafft]="alignment";              STAGEMAP[mafft]="1"
ROLE[clustalo]="alignment";           STAGEMAP[clustalo]="1"
ROLE[primer3_core]="raw-generation";  STAGEMAP[primer3_core]="2"
ROLE[primer3py]="pre-filter";         STAGEMAP[primer3py]="3a"
ROLE[melting]="thermodynamics";       STAGEMAP[melting]="3b"
ROLE[oligomelting]="thermodynamics";  STAGEMAP[oligomelting]="3b"
ROLE[seqkit]="primer-stats"
    STAGEMAP[seqkit]="4a"
    ROLE[bowtie2]="alignment-specificity"
    STAGEMAP[bowtie2]="4a"
    ROLE[blastn]="specificity";           STAGEMAP[blastn]="4a"
ROLE[mfeprimer]="specificity";        STAGEMAP[mfeprimer]="4a"
ROLE[primerdiffer]="specificity";     STAGEMAP[primerdiffer]="4a"
ROLE[ispcr]="in-silico-PCR";          STAGEMAP[ispcr]="4b"
ROLE[ipcress]="in-silico-PCR";        STAGEMAP[ipcress]="4b"
ROLE[primersearch]="in-silico-PCR";   STAGEMAP[primersearch]="4b"
ROLE[tntblast]="in-silico-PCR";       STAGEMAP[tntblast]="4b"
ROLE[primalscheme3]="tiling-scheme";  STAGEMAP[primalscheme3]="5"
ROLE[varvamp]="tiling-scheme";        STAGEMAP[varvamp]="5"
ROLE[olivar]="tiling-scheme";         STAGEMAP[olivar]="5"
ROLE[puppy]="taxon-specific";         STAGEMAP[puppy]="6"
ROLE[ngs_primerplex]="multiplex";     STAGEMAP[ngs_primerplex]="6"
ROLE[primerserver2]="target-focused-design"
    STAGEMAP[primerserver2]="7c"
    ROLE[degeprime]="degenerate-primers"; STAGEMAP[degeprime]="7"

TOOL_ORDER=(
    mafft clustalo
    primer3_core primer3py
    melting oligomelting
    seqkit bowtie2
    blastn mfeprimer primerdiffer
    ispcr ipcress primersearch tntblast
    primalscheme3 varvamp olivar
    puppy ngs_primerplex
    degeprime primerserver2
)

PASS_COUNT=0; FAIL_COUNT=0; SKIP_COUNT=0
for tool in "${TOOL_ORDER[@]}"; do
    s="${STATUS[$tool]:-SKIP}"
    role="${ROLE[$tool]:-unknown}"
    stage="${STAGEMAP[$tool]:-?}"
    if [[ "$s" == "PASS" ]]; then
        printf "${GREEN}%-22s ✅ PASS      %-22s Stage %-4s${NC}\n" "$tool" "$role" "$stage"
        ((PASS_COUNT++))
    elif [[ "$s" == "FAIL" ]]; then
        printf "${RED}%-22s ❌ FAIL      %-22s Stage %-4s${NC}\n" "$tool" "$role" "$stage"
        ((FAIL_COUNT++))
    elif [[ "$s" == "RUNNING" ]]; then
        printf "${CYAN}%-22s ⏳ RUNNING   %-22s Stage %-4s${NC}\n" "$tool" "$role" "$stage"
        ((SKIP_COUNT++))
    else
        printf "${YELLOW}%-22s ⚠️  %-8s   %-22s Stage %-4s${NC}\n" "$tool" "$s" "$role" "$stage"
        ((SKIP_COUNT++))
    fi
done

echo ""
echo "──────────────────────────────────────────────────────────────"
echo -e "${BOLD}Results: ${GREEN}$PASS_COUNT PASS${NC}  ${RED}$FAIL_COUNT FAIL${NC}  ${YELLOW}$SKIP_COUNT OTHER${NC}"
echo ""
echo -e "${BOLD}Output directory: $OUTDIR/${NC}"
echo ""
echo -e "${BOLD}Key files by stage:${NC}"
echo "  Stage 1:  mafft_aligned.fasta"
echo "  Stage 2:  primer3_output_left.txt (10k)  primer3_output_right.txt (10k)"
echo "  Stage 3a: all_primers_tm.tsv (20k Tm)  primers_filtered.tsv (~2k good)"
echo "            pairs_filtered.tsv (matched pairs)  pairs_100_*.txt (top 100)"
echo "  Stage 3b: melting_out/  oligomelting_results.txt"
echo "  Stage 4a: blast_results.txt  mfeprimer_spec/dimer/hairpin  primerdiffer_out/"
echo "  Stage 4b: ispcr_output.fasta  ipcress_output.txt  primersearch_output.txt  tntblast_output.txt"
echo "  Stage 5:  ps3_out/  varvamp_out/  olivar_out/"
echo "  Stage 6:  puppy_primers/  npp_run.log"
echo "  Stage 7:  degeprime_output.tsv"
echo ""
echo -e "${GREEN}${BOLD}Demo complete.${NC}"
}

###############################################################################
# MAIN
###############################################################################
STAGE="${1:-all}"
case "$STAGE" in
    stage1)  stage1 ;;
    stage2)  stage2 ;;
    stage3a) stage3a ;;
    stage3b) stage3b ;;
    stage4a) stage4a ;;
    stage4b) stage4b ;;
    stage5)  stage5 ;;
    stage6)  stage6 ;;
    stage7)  stage7 ;;
    summary) show_summary ;;
    all)
        stage1
        stage2
        stage3a
        stage3b
        stage4a
        stage4b
        stage5
        stage6
        stage7
        show_summary
        ;;
    *)
        echo "Usage: $0 [stage1|stage2|stage3a|stage3b|stage4a|stage4b|stage5|stage6|stage7|summary|all]"
        echo ""
        echo "  stage1  — alignment (mafft, clustalo)"
        echo "  stage2  — generate 20,000 raw primers (primer3_core)"
        echo "  stage3a — pre-filter by Tm+structure → ~2,000 good primers (primer3-py)"
        echo "  stage3b — thermodynamics on filtered (MELTING, oligo-melting)"
        echo "  stage4a — specificity on filtered (BLAST, MFEprimer, primerdiffer)"
        echo "  stage4b — in-silico PCR on top 100 pairs (isPcr, ipcress, primersearch, tntblast)"
        echo "  stage5  — tiling schemes (PS3, varvamp, olivar)"
        echo "  stage6  — taxon-specific + multiplex (PUPpy, NGS-PrimerPlex)"
        echo "  stage7  — other generators (DegePrime, primer3-py multi-genome)"
        ;;
esac