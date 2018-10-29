#!/bin/bash

# This code is executed using parallelization to process multiple samples.

PTHDATA={input.path}
WDIR={output.path}
name={ID}

#####################################################################
### Trim and Filter quality, forward & reverse files (TRIMMOMATIC)
#####################################################################
module load intel/2016 vardictjava/1.0

java -jar {Path}/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $WDIR/$name'R1'.fastq $WDIR/$name'R2'.fastq \
$WDIR/$name'R1_paired.fastq' $WDIR/$name'R1_unpaired.fastq.gz' $WDIR/$name'R2_paired.fastq' $WDIR/$name'R2_unpaired.fastq.gz' \
LEADING:3 TRAILING:3 SLIDINGWINDOW:20:30 MINLEN:220

filesuffix="R1_paired"
filesuffix2="R2_paired"

module purge

####################################
### MOTHUR PIPELINE
####################################
module load gcc/5.3.0 mothur/1.39.5

ffile="$WDIR/$name*$filesuffix"
rfile="$WDIR/$name*$filesuffix2"

mothur "#make.contigs(ffastq=$(echo $ffile.fastq),rfastq=$(echo $rfile.fastq), processors=Autodetect); summary.seqs()"
mothur "#screen.seqs(fasta=$(echo $ffile.trim.contigs.fasta), qfile=$(echo $ffile.trim.contigs.qual), summary=$(echo $ffile.trim.contigs.summary), \
maxambig=0, maxlength=500, maxhomop=10, processors=Autodetect); summary.seqs(); unique.seqs(); count.seqs(); summary.seqs(); \
align.seqs(template=$ALIGN); filter.seqs(); unique.seqs()"

mothur "#pre.cluster(fasta=$(echo $ffile.trim.contigs.good.unique.filter.unique.fasta), \
count=$(echo $ffile.trim.contigs.good.count_table), diffs=1, processors=Autodetect)"

mothur "#chimera.uchime(fasta=$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.fasta), reference=$CHIMERAS, \
dereplicate=t, processors=Autodetect)"
mothur "#remove.seqs(fasta=$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.fasta), \
accnos=$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.ref.uchime.accnos), dups=f)"

mothur "#classify.seqs(fasta=$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.pick.fasta), template=$REFERENCE, \
taxonomy=$TAXONOMY, output=simple, method=wang, processors=Autodetect); \
remove.lineage(taxon=Chloroplast-Mitochondria)"

#########################################
### MOTHUR PIPELINE 2: OTU assignation
#########################################
# Assign sequences to OTUs:
mothur "#cluster.split(fasta=$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.pick.fasta), \
name=$(echo $ffile.trim.contigs.good.unique.filter.names), \
taxonomy=$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.pick.*.wang.taxonomy), \
splitmethod=classify, taxlevel=6, cutoff=0.03, processors=12)"

# Manipulate output .list file from cluster.split to handle it easier (at the end it will show OTU ID and total count for each one): 
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' "$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.pick.opti_mcc.list)" \
> "$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.pick.opti_mcc.transposed.list)"

awk -F '[\s\t,]' '{print $name, NF-1}' "$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.pick.opti_mcc.transposed.list)" \
> "$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.pick.opti_mcc.tr.numbers.list)"

# Find consensus taxonomy for OTUs using .wang.taxonomy file:
mothur "#classify.otu(list=$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.pick.opti_mcc.list), \
count=$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.count_table), \
taxonomy=$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.pick.*.wang.taxonomy), label=0.03)"
sed -i.bak '/unknown;/d' "$(echo $ffile.trim.contigs.good.unique.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy)"

#####################################
### FINAL STEPS
#####################################
##Compress output files in the corresponding folder:
cd $WDIR
module purge
if [ ! -d $WDIR/extract ];then mkdir -p $WDIR/extract;fi
echo "extract files of sample "$name" in directory "$WDIR/extract
tar -cvzf $name.tar.gz $1*

###Remove all non-compressed files after putting them in .tar.gz file. This step can be done manually for more security.
mv $1.tar.gz ../
rm $WDIR/${ID}*
mv ../${ID}.tar.gz .
