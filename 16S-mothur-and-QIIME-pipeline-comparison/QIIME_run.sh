#!/bin/bash

#######################################
##### FUNCTION DESIGNATION: 
#######################################

usage(){
cat << _EUs_
$(basename "$0") [OPTIONS]... [PATHS]... [-p] -- Initiates a pipeline for rRNA 16S microbiome analyses using QIIME.\n
\nDESCRIPTION:\n
	\tWhen using this script, take into account that code is prepared for analysis of 16S rRNA Illumina sequences, as R1 and R2 files in .fastq or\n
	\t.fastq.gz format. Also, modifications on default parameters (trimming conditions, reference database used, chimera removal method...) must be \n
	\tdone INSIDE this code.\n
\nPARAMETERS:\n
	\t-d, --database:\n
		\t\tObligatory parameter. For choosing the database. Valid options are:\n
		\t\t> 'SILVA' - Silva database, release 1.32.\n
		\t\t> 'GREENGENES' - Greengenes database, version aug-2013.\n
	\t-f, --files:\n
		\t\tThis flag indicates that initial files are already present in output dir. If specified, no file transfer or ID file creating will 
		be done.\n
	\t-h, --help:\n
		\t\tPrints this screen.\n
	\t-i, --input:\n
		\t\tPath to the input data. By default is located in home directory.\n
	\t-o, --output:\n
		\t\tPath where output is going to be stored. By default is located in 'HOME/16Sout'.\n
	\t-s, --suffix:\n
		\t\tObligatory parameter. Text string indicating forward & reverse labeling format in sequence file names, singlequoted. \n
		\t\t(e.g., '\.R[1-2]', '_R[1-2]', '_[1-2]'). [Default:'_R[1-2]']\n
	\t-t, --trimming:\n
		\t\tActivates the trimming option (with Trimmomatic-0.36) for cleaning sequences if needed. Default parameters are:\n
		\t\t> Paired-end mode\n
		\t\t> -phred33\n
		\t\t> LEADING:3 TRAILING:3 SLIDINGWINDOW:20:30 MINLEN:220\n
	\t-x, --extension:\n
		\t\tObligatory parameter. Extension of input files, singlequoted (e.g., '.fastq', '.fastq.gz'). [Default:'.fastq']
_EUs_
}

TRIM_fun(){
	module purge
	module load intel/2016 vardictjava/1.0
	mkdir $WDIR/trim_unpaired
	for f in $IDs; do
		java -jar /home/otras/ini/alg/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $WDIR/$f$R1$raw_extns $WDIR/$f$R2$raw_extns \
		$WDIR/$f$R1'_paired'$raw_extns $WDIR/$f$R1'_unpaired'$raw_extns $WDIR/$f$R2'_paired'$raw_extns $WDIR/$f$R2'_unpaired'$raw_extns \
		LEADING:3 TRAILING:3 SLIDINGWINDOW:20:30 MINLEN:220
		mv $WDIR/$f$idsuffix'_unpaired'$raw_extns $WDIR/trim_unpaired
	done
	module purge
}

####################################################################
# VARIABLES:
####################################################################
PTH_data=$HOME
WDIR=$HOME/16Sout
idsuffix='_R[1-2]'
raw_extns='.fastq'
files_cp=
trimming=

UsingID=TotalIDs.txt
export R1=${idsuffix/\[1-2\]/1}
export R2=${idsuffix/\[1-2\]/2}

OPTS=`getopt -o d:f::i:o:s:t::x:h --long database:,files::,input:,output:,suffix:,trimming::,extension:,help -- "$@"`
eval set -- "$OPTS"

while true; do
	case $1 in
		-d | --database)
			export DATABASE=$2; shift 2 ;;
		-f | --files)
			files_cp=1; shift 2 ;;
		-i | --input)
			export PTH_data=$2; shift 2 ;;
		-o | --output)
			export WDIR=$2; shift 2 ;;
		-s | --suffix)
			export idsuffix=$2; shift 2 ;;
		-t | --trimming)
			export trimming=1; shift 2 ;;	
		-x | --extension)
			export raw_extns=$2; shift 2 ;;
		-h | --help)
			echo -e $(usage) | less ; exit ;;
		--) shift ; break ;;
        *) echo "Script definition error! Seems that one or more parameters are not valid..." ; exit 1 ;;
	esac
done

#######################################
##### STARTING SCRIPT: 
#######################################
echo '>>> INPUT DIR: '$PTH_data
echo '>>> OUTPUT DIR: '$WDIR

if [ ! -d $WDIR/Rawdata ]; then mkdir -p $WDIR/Rawdata; fi

### Rawdata shortcuts and ID file creation:
if [ "$files_cp" != "1" ]; then
	ln -s $PTH_data/*$idsuffix$raw_extns $WDIR/Rawdata
	find $WDIR/Rawdata -name *$idsuffix$raw_extns | awk '{ sub(/.*\//,"",$1); sub(/'$idsuffix$raw_extns'/,"",$1); print }' | sort | uniq > $WDIR/TotalIDs.txt
fi
IDs=`awk '{print $1}' $WDIR/TotalIDs.txt`

### -d: SILVA/GREENGENES: 
case $DATABASE in
    SILVA ) echo ">>> DATABASE CHOSEN: SILVA (rel. 132)"
        ;;
    GREENGENES ) echo ">>> DATABASE CHOSEN: GreenGenes (13_8)"
        ;;
    * ) echo "Error: Not a valid database name."
esac

if [ "$DATABASE" = "SILVA" ]; then
	CHIMERADB=$LUSTRE/Silva_DB/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna
	REFERENCEDB=$LUSTRE/Silva_DB/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna
	TAXONOMY=$LUSTRE/Silva_DB/SILVA_132_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_7_levels.txt
elif [ "$DATABASE" = "GREENGENES" ]; then
	CHIMERADB=$LUSTRE/GreenGenes_DB/QIIME_files/gg_13_8_otus/rep_set/99_otus.fasta
	REFERENCEDB=$LUSTRE/GreenGenes_DB/QIIME_files/gg_13_8_otus/rep_set/99_otus.fasta
	TAXONOMY=$LUSTRE/GreenGenes_DB/QIIME_files/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt
fi

lineas=`wc -l $WDIR/TotalIDs.txt | awk -F ' ' '{print $1}'`
echo 'QIIME: Number of animals to analyse = '$lineas

#######################################
##### TRIMMOMATIC: 
#######################################
if [ "$trimming" = "1" ]; then
	echo -e "\n>>> TRIMMING ACTIVATED\n"
	echo $(TRIM_fun)
fi

#######################################
##### QIIME: 
#######################################
echo -e "\n++++++++++++++++QIIME++++++++++++++++\n"
module load qiime

for i in `ls $WDIR/*$idsuffix$raw_extns`; do
	nme=`sed -n $i,1p $WDIR/TotalIDs.txt | tr ' ' '\n' | sed -n 1p`
	# Unzip .fastq.gz data files if they are compressed:
	if ["$raw_extns" = ".fastq.gz" ]; then
		gunzip -c $WDIR/$nme$idsuffix'.fastq.gz' > $WDIR/$nme$idsuffix'.fastq'
	fi
	echo 'Merging sample '$nme
	join_paired_ends.py -j 10 -p 1 -f $WDIR/$nme$R1'_paired.fastq' -r $WDIR/$nme$R2'_paired.fastq' -o $WDIR/QIIME/joined/$nme/ 
	gzip $WDIR/*$nme*
done

## WARNING!! If directory names contain '_', this character must be changed to other character (e.g. '.')
for x in $WDIR/QIIME/joined/*/; do
	nm=`echo $x | xargs -n 1 basename`
	echo $nm
	if [[ $x = *_* ]]; then
		mv $WDIR/QIIME/joined/$nm $WDIR/QIIME/joined/"${nm//_/.}"
	fi
done

## Demultiplex sequences from an input fasta file. Barcodes will be removed from the sequences in the output fasta file by default.
FNAOUT=split_out_prejoin
echo ">> Running Demultiplexing..."

multiple_split_libraries_fastq.py -i $WDIR/QIIME/joined/ -o $WDIR/QIIME/$FNAOUT/ --demultiplexing_method sampleid_by_file \
--include_input_dir_path --remove_filepath_in_name --read_indicator .join
echo -e "Demultiplexing done.\nOutput files ($WDIR/QIIME/$FNAOUT):
- seqs.fna
- histograms.txt
- split_library_log.txt\n"

## Identify chimera with usearch and remove chimer sequences:
export PATH=$PATH:$HOME
ulimit -s unlimited

echo ">> Identifying chimera sequences..."
identify_chimeric_seqs.py -i $WDIR/QIIME/$FNAOUT/seqs.fna -m usearch61 -o $WDIR/QIIME/$FNAOUT/usearch_checked_chimeras/ -r $CHIMERADB
echo "Removing chimera sequences..."
filter_fasta.py -f $WDIR/QIIME/$FNAOUT/seqs.fna -o $WDIR/QIIME/$FNAOUT/seqs_chimeras_filtered.fna \
-s $WDIR/QIIME/$FNAOUT/usearch_checked_chimeras/chimeras.txt -n
echo -e "Chimera filtering done.\nOutput files ($WDIR/QIIME/$FNAOUT):
- seqs_chimeras_filtered.fna\n"

## OTU classifying:
echo ">> Picking OTUs using a closed reference and constructing OTU table..."
BIOMPATH=$WDIR/QIIME/otu_w_tax16S_silva
pick_closed_reference_otus.py -p $HOME/bin/scr/qiime_par_file.txt -i $WDIR/QIIME/$FNAOUT/seqs_chimeras_filtered.fna -o $BIOMPATH/ \
-r $REFERENCEDB -t $TAXONOMY -f
echo -e "Picking OTUs done.\nOutput files ($BIOMPATH):
- otu_table.biom
- uclust_ref_picked_otus/seqs_chimeras_filtered_failures.txt
- uclust_ref_picked_otus/seqs_chimeras_filtered_otus.txt
- uclust_ref_picked_otus/seqs_chimeras_filtered_clusters.uc
- uclust_ref_picked_otus/seqs_chimeras_filtered_otus.log\n"

echo ">> Assigning taxonomy..."
assign_taxonomy.py -i $WDIR/QIIME/$FNAOUT/seqs_chimeras_filtered.fna -o $WDIR/QIIME/ -r $REFERENCEDB -t $TAXONOMY
echo "Assigning taxonomy done."
