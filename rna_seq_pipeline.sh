#!/bin/bash

# ==============================
# INPUTS
# ==============================
SRR_ID=$1
GENOME_FA=$2
ANNOTATION=$3
THREADS=8

# ==============================
# CREATE OUTPUT DIRECTORY
# ==============================
mkdir -p results/${SRR_ID}
cd results/${SRR_ID}

echo "=============================="
echo "🚀 Processing $SRR_ID"
echo "=============================="

# ==============================
# STEP 1: DOWNLOAD FROM ENA (API METHOD)
# ==============================
echo "📥 Fetching download links from ENA..."

FASTQ_URLS=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR_ID}&result=read_run&fields=fastq_ftp" | sed 1d | cut -d',' -f2)

URL1=$(echo $FASTQ_URLS | cut -d';' -f1)
URL2=$(echo $FASTQ_URLS | cut -d';' -f2)

echo "Downloading:"
echo $URL1
echo $URL2

wget -c https://$URL1
wget -c https://$URL2

# Check download
if [ ! -f ${SRR_ID}_1.fastq.gz ] || [ ! -f ${SRR_ID}_2.fastq.gz ]; then
    echo "❌ ERROR: ENA download failed"
    exit 1
fi
# ==============================
# STEP 2: FASTP QC
# ==============================
echo "🧪 Running fastp..."

fastp \
  -i ${SRR_ID}_1.fastq.gz \
  -I ${SRR_ID}_2.fastq.gz \
  -o ${SRR_ID}_1.clean.fastq.gz \
  -O ${SRR_ID}_2.clean.fastq.gz \
  --detect_adapter_for_pe \
  -h ${SRR_ID}_fastp.html \
  -j ${SRR_ID}_fastp.json \
  -w $THREADS

# ==============================
# STEP 3: BUILD INDEX (ONLY ONCE)
# ==============================
cd ../../

if [ ! -f MT_index.00.b.array ]; then
    echo "🧬 Building genome index..."
    subread-buildindex -o /home/ananya/Downloads/TRANS_P/MT_index $GENOME_FA
else
    echo "✅ Index already exists"
fi

cd results/${SRR_ID}

# ==============================
# STEP 4: ALIGNMENT
# ==============================
echo "🔗 Aligning..."

subread-align -T $THREADS \
  -t 1 \
  -i /home/ananya/Downloads/TRANS_P/MT_index \
  -r ${SRR_ID}_1.clean.fastq.gz \
  -R ${SRR_ID}_2.clean.fastq.gz \
  -o ${SRR_ID}.sam

# ==============================
# STEP 5: BAM PROCESSING
# ==============================
echo "📊 Processing BAM..."

samtools view -@ $THREADS -bS ${SRR_ID}.sam > ${SRR_ID}.bam
samtools sort -@ $THREADS ${SRR_ID}.bam -o ${SRR_ID}.sorted.bam
samtools index ${SRR_ID}.sorted.bam

# ==============================
# STEP 5.5: FLAGSTAT QC
# ==============================
echo "📈 Generating alignment stats (flagstat)..."

samtools flagstat ${SRR_ID}.sorted.bam > ${SRR_ID}_flagstat.txt
# ==============================
# STEP 6: FEATURE COUNTS
# ==============================
echo "🧮 Counting reads..."

featureCounts -T $THREADS \
  -p -B \
  -t gene \
  -g locus_tag \
  -a $ANNOTATION \
  -o ${SRR_ID}_counts.txt \
  ${SRR_ID}.sorted.bam


echo "✅ DONE: $SRR_ID"
echo "📁 Output: results/${SRR_ID}/"
