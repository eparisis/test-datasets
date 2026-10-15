#!/usr/bin/env bash
# Rebuild the CyanoSuite mini (-profile test_mini) and mid (-profile test_mid) datasets
# from public ONT reads. See tests/data/README.md for sources and design.
#
#   bash tests/data/make_test_datasets.sh <outdir>
#
# Needs curl and pixi (tools run through `pixi exec`, versions pinned below).
# Already-downloaded reads can be reused with CYANO_READS=/path SPHAE_READS=/path.
# BASE is the URL (or local dir) the files will be served from; it goes into the samplesheets.
set -euo pipefail

OUT=${1:?usage: make_test_datasets.sh <outdir>}
THREADS=${THREADS:-8}
SEED=42
BASE=${BASE:-'${projectDir}/data/test-datasets/'}
mkdir -p "$OUT"/{src,mini,mid}
cd "$OUT"

run() { pixi exec -c conda-forge -c bioconda \
    -s minimap2==2.31 -s samtools==1.24 -s seqkit==2.10.1 -s rasusa==5.1.0 "$@"; }

ENA=https://ftp.sra.ebi.ac.uk/vol1/fastq
NCBI=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi

# --- Sources -------------------------------------------------------------------------
# MicroBench (ENA PRJEB85558), dorado 1.0.0 sup v5.2.0 basecalls
fetch() { [ -s "$2" ] || curl -fsSL "$1" -o "$2"; }
CYANO_READS=${CYANO_READS:-src/ERR15459729.fastq.gz}   # Anabaena (Nostoc) sp. PCC 7120, DSM 107007
SPHAE_READS=${SPHAE_READS:-src/ERR15459737.fastq.gz}   # Sphaerotilus natans DSM 6575
fetch $ENA/ERR154/029/ERR15459729/ERR15459729.fastq.gz "$CYANO_READS"
fetch $ENA/ERR154/037/ERR15459737/ERR15459737.fastq.gz "$SPHAE_READS"

# References: PCC 7120 chromosome + 6 plasmids, S. natans DSM 6575 draft, Nostoc sp. MS1
fetch "$NCBI?db=nuccore&rettype=fasta&id=NC_003272.1,NC_003276.1,NC_003240.1,NC_003267.1,NC_003273.1,NC_003270.1,NC_003241.1" src/pcc7120.fna
fetch https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/689/195/GCF_000689195.1_CLC_Genomics_workbench_5.0/GCF_000689195.1_CLC_Genomics_workbench_5.0_genomic.fna.gz src/sphaerotilus.fna.gz
fetch "$NCBI?db=nuccore&rettype=fasta&id=NZ_AP023441.1,NZ_AP023442.1,NZ_AP023443.1,NZ_AP023444.1,NZ_AP023445.1,NZ_AP023446.1" src/nostoc_ms1.fna

# --- Read-to-reference alignments (primary only, MAPQ >= 20) ---------------------------
align() { # reads ref out.bam
    [ -s "$3" ] || run bash -c "minimap2 -t $THREADS -ax map-ont '$2' '$1' \
        | samtools view -u -F 0x904 -q 20 - | samtools sort -@ $THREADS -o '$3' && samtools index '$3'"
}
align "$CYANO_READS" src/pcc7120.fna src/cyano.bam
align "$SPHAE_READS" src/sphaerotilus.fna.gz src/sphae.bam

# Heterotroph region for mini: the largest S. natans contigs until >= 400 kb (4 contigs, 514 kb)
run seqkit sort -l -r src/sphaerotilus.fna.gz | run seqkit fx2tab -n -i -l \
    | awk '{s+=$2; print $1"\t0\t"$2; if (s>=400000) exit}' > src/sphae_mini.bed

ids() { run samtools view "$@" | cut -f1 | sort -u; }

# --- Mini dataset ----------------------------------------------------------------------
# 1 Mb of the PCC 7120 chromosome (4.5-5.5 Mb, which holds the ribosomal protein cluster, so DAS Tool finds enough single-copy genes) + plasmids delta (55 kb) and zeta (5.6 kb) at ~15x,
# plus 514 kb of S. natans at ~10x (reads >= 2 kb). Different depth and GC (41% vs 70%)
# let the binners separate the two organisms.
ids src/cyano.bam NC_003272.1:4500001-5500000 NC_003273.1 NC_003241.1 > src/mini_cyano.ids
ids -L src/sphae_mini.bed src/sphae.bam > src/mini_sphae.ids
# (rasusa reads its input twice, so it needs a file, not a pipe)
run seqkit grep -f src/mini_cyano.ids "$CYANO_READS" -o src/mini_cyano_all.fastq
run seqkit grep -f src/mini_sphae.ids "$SPHAE_READS" | run seqkit seq -m 2000 -o src/mini_sphae_all.fastq
run rasusa reads -c 15 -g 1061000 -s $SEED -O fastq -Z u src/mini_cyano_all.fastq > src/mini_cyano.fastq
run rasusa reads -c 12 -g 400000 -s $SEED -O fastq -Z u src/mini_sphae_all.fastq > src/mini_sphae.fastq
cat src/mini_cyano.fastq src/mini_sphae.fastq | run seqkit shuffle -s $SEED -o mini/mini_coculture.fastq.gz

# FASTA rows: the same cyano region (reference) and the S. natans fragment (non-cyano fixture)
run seqkit grep -p NC_003272.1 src/pcc7120.fna | run seqkit subseq -r 4500001:5500000 \
    | cat - <(run seqkit grep -p NC_003273.1 -p NC_003241.1 src/pcc7120.fna) \
    | run seqkit seq -i -o mini/mini_cyano.fna.gz
run seqkit grep -f <(cut -f1 src/sphae_mini.bed) src/sphaerotilus.fna.gz \
    | run seqkit seq -i -o mini/mini_noncyano.fna.gz

# --- Mid dataset -----------------------------------------------------------------------
# Whole PCC 7120 (7.21 Mb, chromosome + 6 plasmids) at ~40x mixed with S. natans
# (4.59 Mb) at ~20x, as a non-axenic culture. FASTA row: Nostoc sp. MS1 (RefSeq).
run rasusa reads -c 40 -g 7212000 -s $SEED -O fastq -Z u "$CYANO_READS" > src/mid_cyano.fastq
run rasusa reads -c 20 -g 4593000 -s $SEED -O fastq -Z u "$SPHAE_READS" > src/mid_sphae.fastq
cat src/mid_cyano.fastq src/mid_sphae.fastq | run seqkit shuffle -s $SEED -o mid/mid_coculture.fastq.gz
run seqkit seq -i src/nostoc_ms1.fna -o mid/mid_nostoc_ms1.fna.gz

cat > mini/mini_samplesheet.csv <<EOS
sample,name,fastq,fasta,platform
MINI,PCC7120_coculture,${BASE}mini_coculture.fastq.gz,,ONT
MINI_FA,PCC7120_fragment,,${BASE}mini_cyano.fna.gz,
MINI_FA,Sphaerotilus_fragment,,${BASE}mini_noncyano.fna.gz,
EOS
cat > mid/mid_samplesheet.csv <<EOS
sample,name,fastq,fasta,platform
MID,PCC7120_coculture,${BASE}mid_coculture.fastq.gz,,ONT
MID_FA,Nostoc_MS1,,${BASE}mid_nostoc_ms1.fna.gz,
EOS

rm -f src/*.fastq
run seqkit stats -a -T mini/*.gz mid/*.gz 2>/dev/null | tee stats.tsv
for d in mini mid; do (cd $d && md5sum ./*.gz 2>/dev/null || md5 -r ./*.gz) > $d.md5; done
