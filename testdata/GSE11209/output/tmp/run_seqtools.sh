#!/bin/bash
# This script was generated by BRB-SeqTools 1.0 (Oct 2016). BRB-SeqTools does not take the responsibility for any damage, loss of data resulting from the use of the software.

set -e
export seqtools_samtools_PATH=/opt/SeqTools/bin/samtools-1.3:/opt/SeqTools/bin/samtools-1.3/misc
export PATH=$seqtools_samtools_PATH:$PATH
export seqtools_bowtie_PATH=/opt/SeqTools/bin/bowtie2-2.2.6
export PATH=$seqtools_bowtie_PATH:$PATH
export seqtools_tophat_PATH=/opt/SeqTools/bin/tophat-2.1.0.Linux_x86_64
export PATH=$seqtools_tophat_PATH:$PATH

export seqtools_bcftools_PATH=/opt/SeqTools/bin/bcftools-1.3
export PATH=$seqtools_bcftools_PATH:$PATH

outputDir="/home/brb/testdata/GSE11209/output"
if [ -f "$outputDir"/BugReport.tar.gz ]; then rm "$outputDir"/BugReport.tar.gz; fi;
if ls "$outputDir"/*.count &> /dev/null; then rm "$outputDir"/*.count ; fi
if ls "$outputDir"/*.vcf &> /dev/null; then rm "$outputDir"/*.vcf ; fi
if ls "$outputDir"/*.idx &> /dev/null; then rm "$outputDir"/*.idx ; fi
if [ -d "$outputDir/tmp" ]; then rm -r -f "$outputDir/tmp"; fi

if [ ! -d "$outputDir/tmp" ]; then mkdir "$outputDir/tmp"; fi

if [ -f "$outputDir/tmp/log.txt" ]; then rm "$outputDir/tmp/log.txt"; fi
touch "$outputDir/tmp/log.txt"

echo bowtie: /opt/SeqTools/bin/bowtie2-2.2.6 >> "$outputDir/tmp/log.txt"
echo tophat: /opt/SeqTools/bin/tophat-2.1.0.Linux_x86_64 >> "$outputDir/tmp/log.txt"
echo samtools: /opt/SeqTools/bin/samtools-1.3 >> "$outputDir/tmp/log.txt"
echo HTSeq: `python -c "import HTSeq; print HTSeq.__version__"` >> "$outputDir/tmp/log.txt"

runDGE=true
runVC=true
rerunAligner=true
cleanUp=true

cd "/home/brb/testdata/GSE11209"

# dT_ori 
echo dT_ori `date +'%Y-%m-%d %T'` >> "$outputDir/tmp/log.txt"
echo "dT_ori (1/2)"
if [ ! -f dT_ori/accepted_hits_tophat.bam ] || ([ -f dT_ori/accepted_hits_tophat.bam ] && [ $rerunAligner == "true" ]) ; then
    echo "    running read alignment" | tee -a "$outputDir/tmp/log.txt"
    (tophat2 --no-coverage-search -p 1 -o "dT_ori" \
       -G "/home/brb/testdata/Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Genes/genes.gtf" \
       --transcriptome-index=transcriptome_data/known \
       "/home/brb/testdata/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/Bowtie2Index/genome" \
       SRR002051.fastq ) >> "$outputDir/tmp/log.txt" 2>&1
    if [ -f dT_ori/accepted_hits.bam ]; then
        mv "dT_ori/accepted_hits.bam" "dT_ori/accepted_hits_tophat.bam"
    fi;
fi;

if [ $runDGE == "true" ]; then
    cd "dT_ori"
    echo "    sorting bam file" | tee -a "$outputDir/tmp/log.txt"
    (samtools sort -@ 1 -n accepted_hits_tophat.bam -o sortName.bam) >> "$outputDir/tmp/log.txt" 2>&1
    # samtools view -o sortName.sam sorted.bam
    echo "    creating count file" | tee -a "$outputDir/tmp/log.txt"
    (python -m HTSeq.scripts.count -f bam -s no -a 10 sortName.bam "/home/brb/testdata/Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Genes/genes.gtf" > "/home/brb/testdata/GSE11209/output/dT_ori.count") >> "$outputDir/tmp/log.txt" 2>&1
    if [ $cleanUp == "true" ]; then
        rm -rf `ls | grep -v 'accepted_hits_tophat.bam' | grep -v 'accepted_hits_star.bam' | grep -v 'accepted_hits_bwa.bam' | grep -v 'accepted_hits_usr.bam' | grep -v 'align_summary.txt' | grep -v 'Log.final.out' `
    fi;
    cd ..
fi

if [ $runVC == "true" ]; then
    # SortAndIndexBam. Input: accepted_hits.bam. Output: sorted.bam
    cd "dT_ori"
    echo "    sorting bam file" >> "$outputDir/tmp/log.txt"

    (samtools sort -@ 1 accepted_hits_tophat.bam -o sorted.bam) >> "$outputDir/tmp/log.txt" 2>&1
    echo "    indexing bam file" >> "$outputDir/tmp/log.txt"

    (samtools index sorted.bam) >> "$outputDir/tmp/log.txt" 2>&1
    # MarkDuplication (optional). Input: sort.bam. Output: Libraryname.bam
    echo "    marking duplicates" | tee -a "$outputDir/tmp/log.txt"
    (java -Xmx10g -jar /opt/SeqTools/bin/picard-tools-1.141/picard.jar  MarkDuplicates METRICS_FILE=MarkDudup.metrics INPUT=sorted.bam OUTPUT=dT_ori.bam) >> "$outputDir/tmp/log.txt" 2>&1
    # SamtoolCallVariant. Input: LibraryName.bam. Output: LibraryName_raw.vcf
    echo "    calling variants" | tee -a "$outputDir/tmp/log.txt"
    (samtools mpileup -uf "/home/brb/testdata/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/WholeGenomeFasta/genome.fa" dT_ori.bam | \
          bcftools call -vmO v -o "/home/brb/testdata/GSE11209/output/dT_ori_raw.vcf") >> "$outputDir/tmp/log.txt" 2>&1
    if [ $cleanUp == "true" ]; then
        rm -rf `ls | grep -v 'accepted_hits_tophat.bam' | grep -v 'accepted_hits_star.bam' | grep -v 'accepted_hits_bwa.bam' | grep -v 'accepted_hits_usr.bam'  | grep -v 'align_summary.txt'`
    fi;
    cd ..
    echo "    finish variant call" | tee -a "$outputDir/tmp/log.txt"
    # echo "    dT_ori_raw.vcf" is created | tee -a "$outputDir/tmp/log.txt"
fi

# RH_ori 
echo RH_ori `date +'%Y-%m-%d %T'` >> "$outputDir/tmp/log.txt"
echo "RH_ori (2/2)"
if [ ! -f RH_ori/accepted_hits_tophat.bam ] || ([ -f RH_ori/accepted_hits_tophat.bam ] && [ $rerunAligner == "true" ]) ; then
    echo "    running read alignment" | tee -a "$outputDir/tmp/log.txt"
    (tophat2 --no-coverage-search -p 1 -o "RH_ori" \
       -G "/home/brb/testdata/Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Genes/genes.gtf" \
       --transcriptome-index=transcriptome_data/known \
       "/home/brb/testdata/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/Bowtie2Index/genome" \
       SRR002059.fastq ) >> "$outputDir/tmp/log.txt" 2>&1
    if [ -f RH_ori/accepted_hits.bam ]; then
        mv "RH_ori/accepted_hits.bam" "RH_ori/accepted_hits_tophat.bam"
    fi;
fi;

if [ $runDGE == "true" ]; then
    cd "RH_ori"
    echo "    sorting bam file" | tee -a "$outputDir/tmp/log.txt"
    (samtools sort -@ 1 -n accepted_hits_tophat.bam -o sortName.bam) >> "$outputDir/tmp/log.txt" 2>&1
    # samtools view -o sortName.sam sorted.bam
    echo "    creating count file" | tee -a "$outputDir/tmp/log.txt"
    (python -m HTSeq.scripts.count -f bam -s no -a 10 sortName.bam "/home/brb/testdata/Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Genes/genes.gtf" > "/home/brb/testdata/GSE11209/output/RH_ori.count") >> "$outputDir/tmp/log.txt" 2>&1
    if [ $cleanUp == "true" ]; then
        rm -rf `ls | grep -v 'accepted_hits_tophat.bam' | grep -v 'accepted_hits_star.bam' | grep -v 'accepted_hits_bwa.bam' | grep -v 'accepted_hits_usr.bam' | grep -v 'align_summary.txt' | grep -v 'Log.final.out' `
    fi;
    cd ..
fi

if [ $runVC == "true" ]; then
    # SortAndIndexBam. Input: accepted_hits.bam. Output: sorted.bam
    cd "RH_ori"
    echo "    sorting bam file" >> "$outputDir/tmp/log.txt"

    (samtools sort -@ 1 accepted_hits_tophat.bam -o sorted.bam) >> "$outputDir/tmp/log.txt" 2>&1
    echo "    indexing bam file" >> "$outputDir/tmp/log.txt"

    (samtools index sorted.bam) >> "$outputDir/tmp/log.txt" 2>&1
    # MarkDuplication (optional). Input: sort.bam. Output: Libraryname.bam
    echo "    marking duplicates" | tee -a "$outputDir/tmp/log.txt"
    (java -Xmx10g -jar /opt/SeqTools/bin/picard-tools-1.141/picard.jar  MarkDuplicates METRICS_FILE=MarkDudup.metrics INPUT=sorted.bam OUTPUT=RH_ori.bam) >> "$outputDir/tmp/log.txt" 2>&1
    # SamtoolCallVariant. Input: LibraryName.bam. Output: LibraryName_raw.vcf
    echo "    calling variants" | tee -a "$outputDir/tmp/log.txt"
    (samtools mpileup -uf "/home/brb/testdata/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/WholeGenomeFasta/genome.fa" RH_ori.bam | \
          bcftools call -vmO v -o "/home/brb/testdata/GSE11209/output/RH_ori_raw.vcf") >> "$outputDir/tmp/log.txt" 2>&1
    if [ $cleanUp == "true" ]; then
        rm -rf `ls | grep -v 'accepted_hits_tophat.bam' | grep -v 'accepted_hits_star.bam' | grep -v 'accepted_hits_bwa.bam' | grep -v 'accepted_hits_usr.bam'  | grep -v 'align_summary.txt'`
    fi;
    cd ..
    echo "    finish variant call" | tee -a "$outputDir/tmp/log.txt"
    # echo "    RH_ori_raw.vcf" is created | tee -a "$outputDir/tmp/log.txt"
fi

echo completed `date +'%Y-%m-%d %T'` >> "$outputDir/tmp/log.txt"
