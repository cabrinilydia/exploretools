#!/bin/bash
# ============================================================
# PRIMER DESIGN PIPELINE — FULL SETUP FROM SCRATCH
# Requirements: Linux, internet connection, ~10GB free space
# Time: ~20-30 minutes
# Usage: bash setup_primer_pipeline.sh
# ============================================================

# Always work from the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Do not stop on errors
set +e

echo "========================================"
echo " PRIMER DESIGN PIPELINE SETUP"
echo " Working directory: $SCRIPT_DIR"
echo "========================================"

# ── PREREQUISITE: Miniforge ──────────────────────────────────
# If conda is not installed, install Miniforge first:
#
#   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
#   bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
#   eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
#   conda init bash && source ~/.bashrc
#
# Then re-run this script.

if ! command -v conda &>/dev/null; then
    echo "ERROR: conda not found. Install Miniforge first:"
    echo "  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    echo "  bash Miniforge3-Linux-x86_64.sh -b -p \$HOME/miniforge3"
    echo '  eval "$($HOME/miniforge3/bin/conda shell.bash hook)"'
    echo "  conda init bash && source ~/.bashrc"
    echo "Then re-run this script."
    exit 1
fi

echo "conda found: $(conda --version)"

# ── STEP 1: Conda channels ───────────────────────────────────
echo ""
echo "=== STEP 1: Adding conda channels ==="
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible

# ── STEP 2: Create conda environment ────────────────────────
echo ""
echo "=== STEP 2: Creating conda environment (primer-env) ==="
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

echo "Activating primer-env..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate primer-env
export PATH="$(conda info --base)/envs/primer-env/bin:$PATH"

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

# Upgrade networkx — old version breaks PrimalScheme3 on Python 3.13
pip install "networkx>=3.1" --upgrade

echo "pip packages installed"

# ── STEP 4: Fix PrimerServer2 for Python 3.13 ───────────────
echo ""
echo "=== STEP 4: Patching PrimerServer2 for Python 3.13 ==="

SITE_PKGS=$(python3 -c 'import site; print(site.getsitepackages()[0])')

# Fix 4a: distutils removed in Python 3.12+
sed -i 's/from distutils.version import LooseVersion/from packaging.version import Version as LooseVersion/' \
    "$SITE_PKGS/primerserver2/cmd/primertool.py"
echo "  primertool.py patched"

# Fix 4b: primer3.bindings API changed in primer3-py v2
python3 - << 'INNEREOF'
import re, site
filepath = site.getsitepackages()[0] + "/primerserver2/core/design_primer.py"
with open(filepath) as f:
    content = f.read()
content = content.replace('primer3.bindings.setP3Globals(', '_p3_global = ')
content = content.replace('primer3.bindings.designPrimers(', 'primer3.design_primers(')
content = re.sub(
    r'(primer3\.design_primers\(\s*\{[^}]+\}\s*)\)',
    r'\1, _p3_global)',
    content,
    flags=re.DOTALL
)
# Fix stray ) left by replacement
content = content.replace('_p3_global = p3_settings)\n', '_p3_global = p3_settings\n')
with open(filepath, 'w') as f:
    f.write(content)
print("  design_primer.py patched")
INNEREOF

# Fix 4c: PRIMER_PAIR_NUM_RETURNED_FINAL key removed in new API
sed -i 's/PRIMER_PAIR_NUM_RETURNED_FINAL/PRIMER_PAIR_NUM_RETURNED/g' \
    "$SITE_PKGS/primerserver2/core/output.py"

# Fix 4d: PRIMER_PAIR_AMPLICON_NUM_RANK key missing in design mode
python3 - << 'INNEREOF'
import site
filepath = site.getsitepackages()[0] + "/primerserver2/core/output.py"
with open(filepath) as f:
    content = f.read()
old = "raw_rank = primers[f'PRIMER_PAIR_AMPLICON_NUM_RANK_{amplicon_rank}']"
new = "raw_rank = primers.get(f'PRIMER_PAIR_AMPLICON_NUM_RANK_{amplicon_rank}', amplicon_rank)"
content = content.replace(old, new)
with open(filepath, 'w') as f:
    f.write(content)
print("  output.py patched")
INNEREOF

echo "PrimerServer2 fully patched"

# ── STEP 5: Download MFEprimer binary ───────────────────────
echo ""
echo "=== STEP 5: Downloading MFEprimer binary ==="
mkdir -p ~/bin

wget -q "https://github.com/quwubin/MFEprimer-3.0/releases/download/v4.2.4/mfeprimer-4.2.4-linux-amd64.gz" \
    -O ~/bin/mfeprimer.gz && \
    gzip -df ~/bin/mfeprimer.gz && \
    chmod +x ~/bin/mfeprimer && \
    echo "  MFEprimer downloaded" || \
    echo "  WARNING: MFEprimer download failed — check URL manually"

export PATH="$HOME/bin:$PATH"
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc

# ── STEP 6: Clone git tools ──────────────────────────────────
echo ""
echo "=== STEP 6: Cloning git tools ==="

# DegePrime
if [[ -d ~/DegePrime ]]; then
    echo "  DegePrime already exists, skipping"
else
    git clone https://github.com/EnvGen/DegePrime ~/DegePrime
    echo "  DegePrime cloned"
fi

# NGS-PrimerPlex
if [[ -d ~/NGS-PrimerPlex ]]; then
    echo "  NGS-PrimerPlex already exists, skipping"
else
    git clone https://github.com/aakechin/NGS-PrimerPlex ~/NGS-PrimerPlex
    echo "  NGS-PrimerPlex cloned"
fi
TMPDIR=~/tmp pip install -r ~/NGS-PrimerPlex/linux_requirements.txt -q
echo "  NGS-PrimerPlex dependencies installed"

# PUPpy
if [[ -d ~/PUPpy ]]; then
    echo "  PUPpy already exists, skipping"
else
    git clone https://github.com/Dreycey/PUPpy ~/PUPpy
    echo "  PUPpy cloned"
fi
echo 'export PATH="$HOME/PUPpy/scripts:$PATH"' >> ~/.bashrc
export PATH="$HOME/PUPpy/scripts:$PATH"

# ── STEP 7: MELTING (Java) ───────────────────────────────────
echo ""
echo "=== STEP 7: Downloading MELTING (Java Tm calculator) ==="
if [[ -d ~/MELTING5.2.0 ]]; then
    echo "  MELTING already exists, skipping"
else
    wget -q "https://sourceforge.net/projects/melting/files/melting5/MELTING5.2.0.zip" \
        -O /tmp/MELTING5.2.0.zip && \
        unzip -q /tmp/MELTING5.2.0.zip -d ~/ && \
        rm /tmp/MELTING5.2.0.zip && \
        echo "  MELTING installed" || \
        echo "  WARNING: MELTING download failed"
fi

# ── STEP 8: Create folder structure ─────────────────────────
echo ""
echo "=== STEP 8: Creating folder structure ==="
mkdir -p "$SCRIPT_DIR/tool_outputs_final"
mkdir -p "$SCRIPT_DIR/bacterial_cds/target"
mkdir -p "$SCRIPT_DIR/bacterial_cds/nontarget"
echo "  Folders created in $SCRIPT_DIR"

# ── STEP 9: Download test data from NCBI ────────────────────
echo ""
echo "=== STEP 9: Downloading 5 SARS-CoV-2 genomes from NCBI ==="

if [[ -f "$SCRIPT_DIR/test5.fasta" ]]; then
    COUNT=$(grep -c ">" "$SCRIPT_DIR/test5.fasta")
    echo "  test5.fasta already exists ($COUNT genomes), skipping"
else
    efetch -db nucleotide \
        -id MN908947.3,MN985325.1,MN988713.1,MN938384.1,MN975262.1 \
        -format fasta > "$SCRIPT_DIR/test5.fasta"
    COUNT=$(grep -c ">" "$SCRIPT_DIR/test5.fasta")
    echo "  Downloaded $COUNT genomes → test5.fasta"
fi

# ── STEP 10: Download bacterial CDS for PUPpy ───────────────
echo ""
echo "=== STEP 10: Downloading bacterial CDS (E. coli + Salmonella) ==="

if [[ -s "$SCRIPT_DIR/bacterial_cds/target/EcoliK12_cds.fna" ]]; then
    echo "  EcoliK12_cds.fna already exists, skipping"
else
    efetch -db nuccore \
        -query "Escherichia coli K-12 MG1655[organism]" \
        -format fasta_cds_na 2>/dev/null | head -500 > /tmp/ecoli_raw.fna
    python3 -c "
from Bio import SeqIO
seqs = list(SeqIO.parse('/tmp/ecoli_raw.fna', 'fasta'))[:20]
SeqIO.write(seqs, '$SCRIPT_DIR/bacterial_cds/target/EcoliK12_cds.fna', 'fasta')
print(f'  E. coli: {len(seqs)} CDS sequences saved')
"
fi

if [[ -s "$SCRIPT_DIR/bacterial_cds/nontarget/SalmonellaTM_cds.fna" ]]; then
    echo "  SalmonellaTM_cds.fna already exists, skipping"
else
    efetch -db nuccore \
        -query "Salmonella enterica Typhimurium[organism]" \
        -format fasta_cds_na 2>/dev/null | head -500 > /tmp/salmonella_raw.fna
    python3 -c "
from Bio import SeqIO
seqs = list(SeqIO.parse('/tmp/salmonella_raw.fna', 'fasta'))[:20]
SeqIO.write(seqs, '$SCRIPT_DIR/bacterial_cds/nontarget/SalmonellaTM_cds.fna', 'fasta')
print(f'  Salmonella: {len(seqs)} CDS sequences saved')
"
fi

# ── STEP 11: Verify all tools ────────────────────────────────
echo ""
echo "=== STEP 11: Verifying tools ==="
echo ""

check_tool() {
    if command -v "$1" &>/dev/null; then
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

if [[ -x "$HOME/bin/mfeprimer" ]]; then
    echo "  ✅ mfeprimer"
else
    echo "  ❌ mfeprimer NOT FOUND"
fi

echo ""
python3 -c "import primer3;     print('  ✅ primer3-py')"   2>/dev/null || echo "  ❌ primer3-py"
python3 -c "import dask;        print('  ✅ dask')"         2>/dev/null || echo "  ❌ dask"
python3 -c "import pyarrow;     print('  ✅ pyarrow')"      2>/dev/null || echo "  ❌ pyarrow"
python3 -c "import pysam;       print('  ✅ pysam')"        2>/dev/null || echo "  ❌ pysam"
python3 -c "import seaborn;     print('  ✅ seaborn')"      2>/dev/null || echo "  ❌ seaborn"
python3 -c "import networkx; print('  ✅ networkx', networkx.__version__)" 2>/dev/null || echo "  ❌ networkx"

echo ""
echo "========================================"
echo " SETUP COMPLETE"
echo "========================================"
echo ""
echo "Next steps:"
echo "  1. conda activate primer-env"
echo "  2. cd $SCRIPT_DIR"
echo "  3. bash exploretoolsfinal.sh all 2>&1 | tee demo_run_final.log"
