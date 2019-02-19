# BAM Directory
BD=06_PMC_markers_bam
SUFFIX=bam
# Stats directory
SD=06_PMC_markers_bam_stats

mkdir $SD

for F in ${BD}/*.${SUFFIX}
do
	NODIR=$(basename ${F})
	BASE=`echo "$NODIR" | cut -d'.' -f1`
	#run samtools flagstat
	if [ ! -f  ${SD}/${BASE}.stats ]; then
		samtools flagstat $F > ${SD}/${BASE}.stats
	fi
done

######################################

# BAM Directory
BD=06_PMC_markers_group_bam
SUFFIX=bam
# Stats directory
SD=06_PMC_markers_group_bam_stats

mkdir $SD

for F in ${BD}/*.${SUFFIX}
do
	NODIR=$(basename ${F})
	BASE=`echo "$NODIR" | cut -d'.' -f1`
	#run samtools flagstat
	if [ ! -f  ${SD}/${BASE}.stats ]; then
		samtools flagstat $F > ${SD}/${BASE}.stats
	fi
done

######################################

# BAM Directory
BD=06_PMC_markers_subclass_bam
SUFFIX=bam
# Stats directory
SD=06_PMC_markers_subclass_bam_stats

mkdir $SD

for F in ${BD}/*.${SUFFIX}
do
	NODIR=$(basename ${F})
	BASE=`echo "$NODIR" | cut -d'.' -f1`
	#run samtools flagstat
	if [ ! -f  ${SD}/${BASE}.stats ]; then
		samtools flagstat $F > ${SD}/${BASE}.stats
	fi
done
