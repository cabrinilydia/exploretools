#!/bin/bash
# Get directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# ============================================================
# PRIMER DESIGN PIPELINE — FULL SETUP FROM SCRATCH
# Run on fresh Linux machine with Miniforge installed
# Usage: bash setup_primer_pipeline.sh
# Time: ~20-30 minutes
# ============================================================

set +e

echo "========================================"
echo " PRIMER DESIGN PIPELINE SETUP"
echo "========================================"

# ── STEP 1: Conda channels ───────────────────────────────────
echo ""
echo "=== STEP 1: Adding conda channels ==="
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible

# ── STEP 2: Create conda environment ────────────────────────
echo ""
echo "=== STEP 2: Creating conda environment ==="
conda create -n primer-env python=3.13 \
    biopython \
    mafft \
    clustalo \
    primer3 \
    blast \
    bowtie2 \
    samtools \
    seqkit \
    emboss \
    exonerate \
    ispcr \
    tntblast \
    olivar \
    openjdk \
    parasail-python \
    bedops \
    ninja \
    numpy \
    entrez-direct \
    mmseqs2 -y

echo "Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate primer-env

# ── STEP 3: pip packages ─────────────────────────────────────
echo ""
echo "=== STEP 3: Installing pip packages ==="
mkdir -p ~/tmp

TMPDIR=~/tmp pip install \
    primer3-py \
    primalscheme3 \
    varvamp \
    oligo-melting \
    primerdiffer \
    packaging \
    primerserver2 \
    pysam \
    colorama \
    dask \
    distributed \
    pyarrow \
    seaborn \
    matplotlib \
    scipy

pip install "networkx>=3.1" --upgrade 2>/dev/null

# ── STEP 4: Fix PrimerServer2 Python 3.13 ───────────────────
echo ""
echo "=== STEP 4: Fixing PrimerServer2 for Python 3.13 ==="
sed -i 's/from distutils.version import LooseVersion/from packaging.version import Version as LooseVersion/' \
    $(python3 -c 'import site; print(site.getsitepackages()[0])')/primerserver2/cmd/primertool.py
echo "PrimerServer2 patched"

# Patch design_primer.py for primer3-py v2 API
python3 - << INNEREOF
import re
filepath = "$(python3 -c 'import site; print(site.getsitepackages()[0])')/primerserver2/core/design_primer.py"
with open(filepath) as f:
    content = f.read()
content = content.replace('primer3.bindings.setP3Globals(', '_p3_global = ')
content = content.replace('primer3.bindings.designPrimers(', 'primer3.design_primers(')
content = re.sub(r'(primer3\.design_primers\(\s*\{[^}]+\}\s*)\)', r'\1, _p3_global)', content, flags=re.DOTALL)
with open(filepath, 'w') as f:
    f.write(content)
print("design_primer.py patched")
INNEREOF

# Patch output.py
sed -i 's/PRIMER_PAIR_NUM_RETURNED_FINAL/PRIMER_PAIR_NUM_RETURNED/g' \
    $(python3 -c 'import site; print(site.getsitepackages()[0])')/primerserver2/core/output.py
echo "output.py patched"

# ── STEP 5: Download binaries ────────────────────────────────
mkdir -p ~/bin

wget -q https://github.com/quwubin/MFEprimer-3.0/releases/download/v4.2.4/mfeprimer-4.2.4-linux-amd64.gz \
  -O ~/bin/mfeprimer.gz

rm -f ~/bin/mfeprimer
gzip -df ~/bin/mfeprimer.gz
chmod +x ~/bin/mfeprimer

export PATH="$HOME/bin:$PATH"
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc

rsync -a hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/isPcr ~/bin/ 2>/dev/null || \
    echo "WARNING: isPcr download failed — retry manually"
chmod +x ~/bin/isPcr 2>/dev/null || true
echo "isPcr downloaded"

echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc

# ── STEP 6: Git clone tools ──────────────────────────────────
echo ""
echo "=== STEP 6: Cloning git tools ==="

# DegePrime
if [[ -d ~/DegePrime ]]; then
    echo "DegePrime already exists, skipping"
else
    git clone https://github.com/EnvGen/DegePrime ~/DegePrime
    echo "DegePrime cloned"
fi

# NGS-PrimerPlex
if [[ -d ~/NGS-PrimerPlex ]]; then
    echo "NGS-PrimerPlex already exists, skipping"
else
    git clone https://github.com/aakechin/NGS-PrimerPlex ~/NGS-PrimerPlex
    echo "NGS-PrimerPlex cloned"
fi
TMPDIR=~/tmp pip install -r ~/NGS-PrimerPlex/linux_requirements.txt -q

# PUPpy
if [[ -d ~/PUPpy ]]; then
    echo "PUPpy already exists, skipping"
else
    git clone https://github.com/Dreycey/PUPpy ~/PUPpy
    echo "PUPpy cloned"
fi
echo 'export PATH="$HOME/PUPpy/scripts:$PATH"' >> ~/.bashrc

# ── STEP 7: MELTING ──────────────────────────────────────────
echo ""
echo "=== STEP 7: Downloading MELTING (Java) ==="
if [[ -d ~/MELTING5.2.0 ]]; then
    echo "MELTING already exists, skipping"
else
    wget -q https://sourceforge.net/projects/melting/files/melting5/MELTING5.2.0.zip
    unzip -q MELTING5.2.0.zip -d ~/
    rm MELTING5.2.0.zip
    echo "MELTING installed"
fi

# ── STEP 8: Working directory ────────────────────────────────
echo ""
echo "=== STEP 8: Creating working directory ==="
mkdir -p ~/primer/primer-framework
cd ~/primer/primer-framework
mkdir -p tool_outputs_final
mkdir -p bacterial_cds/target
mkdir -p bacterial_cds/nontarget
echo "Working directory: ~/primer/primer-framework"

# ── STEP 9: Download SARS-CoV-2 test data from NCBI ─────────
echo ""
echo "=== STEP 9: Downloading 5 SARS-CoV-2 genomes from NCBI ==="
efetch -db nucleotide \
    -id MN908947.3,MN985325.1,MN988713.1,MN938384.1,MN975262.1 \
    -format fasta > test5.fasta

COUNT=$(grep -c ">" test5.fasta)
echo "Downloaded $COUNT genomes"
if [[ "$COUNT" -ne 5 ]]; then
    echo "WARNING: Expected 5 genomes, got $COUNT"
fi

# ── STEP 10: Download bacterial CDS for PUPpy ───────────────
echo ""
echo "=== STEP 10: Downloading bacterial CDS (E. coli + Salmonella) ==="

# Download and keep only first 20 sequences (small subset for demo)
# Get script directory so files go to the right place
BACT_TARGET="$SCRIPT_DIR/bacterial_cds/target"
BACT_NONTARGET="$SCRIPT_DIR/bacterial_cds/nontarget"
mkdir -p "$BACT_TARGET" "$BACT_NONTARGET"

efetch -db nuccore \
    -query "Escherichia coli K-12 MG1655[organism]" \
    -format fasta_cds_na 2>/dev/null | head -500 > /tmp/ecoli_raw.fna

python3 -c "
from Bio import SeqIO
seqs = list(SeqIO.parse('/tmp/ecoli_raw.fna', 'fasta'))[:20]
SeqIO.write(seqs, '$BACT_TARGET/EcoliK12_cds.fna', 'fasta')
print(f'E. coli: {len(seqs)} CDS sequences saved')
"

efetch -db nuccore \
    -query "Salmonella enterica Typhimurium[organism]" \
    -format fasta_cds_na 2>/dev/null | head -500 > /tmp/salmonella_raw.fna

python3 -c "
from Bio import SeqIO
seqs = list(SeqIO.parse('/tmp/salmonella_raw.fna', 'fasta'))[:20]
SeqIO.write(seqs, '$BACT_NONTARGET/SalmonellaTM_cds.fna', 'fasta')
print(f'Salmonella: {len(seqs)} CDS sequences saved')
"

# ── STEP 11: Reload PATH ─────────────────────────────────────
echo ""
echo "=== STEP 11: Reloading PATH ==="
source ~/.bashrc
export PATH=$HOME/bin:$HOME/PUPpy/scripts:$PATH

# ── STEP 12: Verify all tools ────────────────────────────────
echo ""
echo "=== STEP 12: Verifying tools ==="
echo ""

check_tool() {
    if command -v $1 &>/dev/null; then
        echo "  ✅ $1"
    else
        echo "  ❌ $1 NOT FOUND"
    fi
}

check_tool mafft
check_tool clustalo
check_tool primer3_core
check_tool blastn
check_tool bowtie2
check_tool samtools
check_tool seqkit
check_tool isPcr
check_tool tntblast
check_tool primersearch
check_tool ipcress
check_tool oligomelting
check_tool varvamp
check_tool primerdesign.py
check_tool primalscheme3
check_tool primertool
check_tool puppy-align
check_tool puppy-primers
check_tool ~/bin/mfeprimer

echo ""
python3 -c "import primer3; print('  ✅ primer3-py')" 2>/dev/null || echo "  ❌ primer3-py"
python3 -c "import dask; print('  ✅ dask')" 2>/dev/null || echo "  ❌ dask"
python3 -c "import pyarrow; print('  ✅ pyarrow')" 2>/dev/null || echo "  ❌ pyarrow"
python3 -c "import pysam; print('  ✅ pysam')" 2>/dev/null || echo "  ❌ pysam"
python3 -c "import seaborn; print('  ✅ seaborn')" 2>/dev/null || echo "  ❌ seaborn"

echo ""
echo "========================================"
echo " SETUP COMPLETE"
echo "========================================"
echo ""
echo "Next steps:"
echo "  1. Put exploretoolsfinal.sh in the SAME folder as this script: $SCRIPT_DIR"
echo "  2. conda activate primer-env"
echo "  3. cd $SCRIPT_DIR"
echo "  4. bash exploretoolsfinal.sh all"
echo ""
# Fix output.py for PrimerServer2
python3 - << INNEREOF
filepath = "$(python3 -c 'import site; print(site.getsitepackages()[0])')/primerserver2/core/output.py"
with open(filepath) as f:
    content = f.read()
old = "raw_rank = primers[f'PRIMER_PAIR_AMPLICON_NUM_RANK_{amplicon_rank}']"
new = "raw_rank = primers.get(f'PRIMER_PAIR_AMPLICON_NUM_RANK_{amplicon_rank}', amplicon_rank)"
content = content.replace(old, new)
with open(filepath, 'w') as f:
    f.write(content)
print("output.py fixed")
INNEREOF
