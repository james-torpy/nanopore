#!/bin/bash

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/single_cell/nanopore"
in_dir="$project_dir/FM/subset"
temp_dir="$in_dir/temp"

READS=$( ls $in_dir/*fastq | sort -h -t "_" -k 3 | sed -n ${SGE_TASK_ID}p  )
T_READS=${READS##*/}
THREADS=4

#module load marsmi/gcc-4.8.4
module load evaben/gcc/gcc-7.3.0/7.3.0 
module load marsmi/nanopore/canu 
module load marsmi/nanopore/minimap2/2.7-r667-dirty
module load marsmi/nanopore/racon

if [[ ! -d logs ]] ; then mkdir logs ;  fi  
#if [[ ! -d bams ]] ; then mkdir bams ;  fi
if [[ ! -d fasta ]] ; then mkdir fasta ; fi
if [[ ! -d temp ]] ; then mkdir temp ; fi

>&2 echo "[*] Launching contig assembly for "$READS
canu  -d $temp_dir/canu \
  -p test \
  -useGrid=false \
  -Overlapper=minimap \
  -batMemory=56 \
  -minReadLength=500 \
  -minOverlapLength=100 \
  -genomeSize=1k \
  -cnsConcurrency=2 \
  -cnsErrorRate=0.25 \
  -maxThreads=${THREADS} \
  -minThreads=${THREADS} \
  -nanopore-raw ${READS} 2> logs/canu.${SGE_TASK_ID}.log \
  -stopOnReadQuality=false

DRAFT=$temp_dir/canu/test.contigs.fasta
if [[ ! -s ${DRAFT} ]] ; then 
  echo "[*] No contigs assembled... exiting" 
  qstat -j ${JOB_ID} | grep usage
  exit 0 
fi

cp $DRAFT fasta/${T_READS%*.fastq}_contigs_canu.${SGE_TASK_ID}.fa
echo "[*] Launching consensus calculation"
for ROUND in {01..04}; do
  echo "Running round ${ROUND} consensus..."
  READS2TIGS=$temp_dir/reads2contigs_${ROUND}.paf
  NEWDRAFT=$temp_dir/racon_${ROUND}.fasta
  minimap2 -k 15 -t ${THREADS} ${DRAFT} ${READS} > ${READS2TIGS}
  racon -w 200 -m 8 -x -1 -g -4  -t ${THREADS} -q -1 ${READS} ${READS2TIGS} ${DRAFT} > ${NEWDRAFT}
  #rm ${DRAFT}
  DRAFT=${NEWDRAFT}
done 2> logs/consensus.${SGE_TASK_ID}.log
cp $DRAFT fasta/${T_READS%*.fastq}_contigs_racon.${SGE_TASK_ID}.fa

# echo "[*] mapping reads to consensus assembly"
# BAM=$(pwd)/bams/${T_READS%*.fastq}_map2contig.${SGE_TASK_ID}.bam
# ( minimap2 -k 15 -a -t ${THREADS}  ${DRAFT} ${READS} | samtools view -b - | samtools sort -@ 3 -o ${BAM} ) 2>&1 >> logs/minimap.${SGE_TASK_ID}.log
# samtools index ${BAM} 2>&1 > logs/samtools.${SGE_TASK_ID}.log
mv ${DRAFT} fasta/${T_READS%*.fastq}_contigs_racon.${SGE_TASK_ID}.fa
rm $temp_dir/racon* $temp_dir/reads2contigs*

echo "All done! Usage: "
qstat -j ${JOB_ID} | grep usage#