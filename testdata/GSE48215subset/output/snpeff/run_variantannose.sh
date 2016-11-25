#!/bin/bash
seqtools_snpeff=/opt/SeqTools/bin/snpEff
export bdge_bcftools_PATH=/opt/SeqTools/bin/bcftools-1.3
export PATH=$bdge_bcftools_PATH:$PATH

export bdge_htslib_PATH=/opt/SeqTools/bin/samtools-1.3/htslib-1.3
export PATH=$bdge_htslib_PATH:$PATH

RcodeDir="$(dirname /home/brb/seqtools-dl/SeqTools)/code"
if [ ! -d "$RcodeDir" ]; then exit 8; fi;inputVCFDir="/home/brb/testdata/GSE48215subset/output"
cosmicVCF="/home/brb/variantAnnoDatabase/cosmic/GRCh37/CosmicCodingMuts.vcf.gz"
dbSNPVCF="/home/brb/variantAnnoDatabase/dbsnp/GRCh37/common_all_20160601.vcf.gz"
genomeRef="/home/brb/testdata/hg19/chr1.fa"
minQual=20
minReadDepth=5
minMapQual=1
genomeVer=GRCh37
outputDir="/home/brb/testdata/GSE48215subset/output/snpeff"
if [ "$genomeVer" == "GRCh37" ]; then
genomeVer="GRCh37.75"
fi
if [ "$genomeVer" == "GRCh38" ]; then
genomeVer="GRCh38.82"
fi
if [ ! -d "$outputDir" ]; then mkdir "$outputDir"; fi
if ls "$outputDir"/*.txt >/dev/null; then rm "$outputDir"/*.txt ; fi
if ls "$outputDir"/*_annotated.vcf >/dev/null; then rm "$outputDir"/*_annotated.vcf ; fi
if ls "$outputDir"/*.pdf >/dev/null; then rm "$outputDir"/*.pdf ; fi
if ls "$outputDir"/*.html >/dev/null; then rm "$outputDir"/*.html ; fi
if ls "$outputDir"/run_variantanno.sh >/dev/null; then rm "$outputDir"/run_variantanno.sh ; fi
if [ -d "$outputDir/tmp" ]; then rm -r -f "$outputDir/tmp";  fi
mkdir "$outputDir/tmp"
if [ -f "$outputDir/tmp/log.txt" ]; then rm "$outputDir/tmp/log.txt"; fi
echo started `date +'%Y-%m-%d %T'` >> "$outputDir/tmp/log.txt"
#Automatically download academic version of dbNSFP database files from BRB-SeqTools googledrive links
userDir="/home/brb"
chmod a+x "$RcodeDir/gdown.pl"
if [ ! -d "$userDir/variantAnnoDatabase" ]; then mkdir "$userDir/variantAnnoDatabase";  fi
if [ ! -d "$userDir/variantAnnoDatabase/dbNSFP" ]; then mkdir "$userDir/variantAnnoDatabase/dbNSFP";  fi
if [ "$genomeVer" == "hg19" ] || [ "$genomeVer" == "GRCh37.75" ]; then
file1="$userDir/variantAnnoDatabase/dbNSFP/dbNSFP3.2a_hg19.txt.gz.tbi"
   if [ ! -f "$file1" ]; then
      ("$RcodeDir/./gdown.pl" 'https://drive.google.com/uc?export=download&id=0B4SLNoZhyUgoMEt2b1N2X3c4TmM' "$file1") 2>&1 | tee -a "$outputDir/tmp/log.txt"
   fi
   file2="$userDir/variantAnnoDatabase/dbNSFP/dbNSFP3.2a_hg19.txt.gz"
   if [ ! -f "$file2" ]; then
      ("$RcodeDir/./gdown.pl" 'https://drive.google.com/uc?export=download&id=0B4SLNoZhyUgod3ZmT1pUM21VSEE' "$file2") 2>&1 | tee -a "$outputDir/tmp/log.txt"
   fi
fi
if [ "$genomeVer" == "hg38" ] || [ "$genomeVer" == "GRCh38.82" ]; then
   file1="$userDir/variantAnnoDatabase/dbNSFP/dbNSFP3.2a_hg38_sorted.txt.gz.tbi"
    if [ ! -f "$file1" ]; then
     ("$RcodeDir/./gdown.pl" 'https://drive.google.com/uc?export=download&id=0BxfUYar-G1WuRGVIY05JSjVlaFU' "$file1") 2>&1 | tee -a "$outputDir/tmp/log.txt"
   fi
   file2="$userDir/variantAnnoDatabase/dbNSFP/dbNSFP3.2a_hg38_sorted.txt.gz"
   if [ ! -f "$file2" ]; then
     ("$RcodeDir/./gdown.pl" 'https://drive.google.com/uc?export=download&id=0BxfUYar-G1WueXp6aEhvU1dwbTQ' "$file2") 2>&1 | tee -a "$outputDir/tmp/log.txt"
   fi
fi
dbnsfpFile="$file2"
count2=$(head -1 "$genomeRef" | grep "chr" | wc -l)
if [[ $dbSNPVCF == *.vcf ]]; then count3=$(grep -v '^#' "$dbSNPVCF" | head -1 | grep "chr" | wc -l); fi;
if [[ $dbSNPVCF == *.vcf.gz ]]; then 
    count3=$(zgrep "^chr[0-9]" "$dbSNPVCF" | head -1 | wc -l)
fi;
for inputVCF in $inputVCFDir/*.vcf
do
tmpname="$(basename $inputVCF)"
outputFileNameShare="${tmpname::-4}"
tmpfd="tmp_$outputFileNameShare"
if [ ! -d "$outputDir/tmp/$tmpfd" ]; then mkdir "$outputDir/tmp/$tmpfd";  fi

echo "running sample: $inputVCF" >> "$outputDir/tmp/log.txt"
if [ -f "$outputDir/tmp/$tmpfd/log.txt" ]; then rm "$outputDir/tmp/log.txt"; fi
bcftools filter -i"QUAL >= $minQual && DP >= $minReadDepth && MQ >= $minMapQual" "$inputVCF" > "$outputDir/tmp/$tmpfd/filtered.vcf"
read num0a <<< $(grep -v '#' "$outputDir/tmp/$tmpfd/filtered.vcf" | wc -l)
#Added an error handler if no varaints remain after one raw VCF has been filtered by the DP, MAQ, QUAL criteria.
if (( $num0a == 0 ))
then 
exit 3
fi
(bcftools norm -m-both -o "$outputDir/tmp/$tmpfd/splitted.vcf" "$outputDir/tmp/$tmpfd/filtered.vcf") 2>&1 | tee -a "$outputDir/tmp/log.txt"
count=$(grep -v '#' "$inputVCF" | cut -f1 | head -1 | grep 'chr' | wc -l)
if [ $count == "1" ] && [ $count2 == "0" ]
then
    perl -pe 's/^chr//' "$outputDir/tmp/$tmpfd/splitted.vcf" > "$outputDir/tmp/$tmpfd/splitted.nochr.vcf"
    (bcftools norm -f $genomeRef -o "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" "$outputDir/tmp/$tmpfd/splitted.nochr.vcf") 2>&1 | tee -a "$outputDir/tmp/log.txt"
    if [ $count3 == "1" ] 
    then
        perl -pe 's/^([^#])/chr\1/' "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" > "$outputDir/tmp/$tmpfd/leftnormalized.vcf"
    else
    mv "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" "$outputDir/tmp/$tmpfd/leftnormalized.vcf"
    fi
fi 
if [ $count == "0" ] && [ $count2 == "1" ]
then
    perl -pe 's/^([^#])/chr\1/' "$outputDir/tmp/$tmpfd/splitted.vcf" > "$outputDir/tmp/$tmpfd/splitted.addchr.vcf"
    grep -v "^chrMT" "$outputDir/tmp/$tmpfd/splitted.addchr.vcf" > "$outputDir/tmp/$tmpfd/splitted.addchr2.vcf"
    (bcftools norm -f $genomeRef -o "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" "$outputDir/tmp/$tmpfd/splitted.addchr2.vcf") 2>&1 | tee -a "$outputDir/tmp/log.txt" 
    if [ $count3 == "0" ]
    then
    grep "^MT" "$outputDir/tmp/$tmpfd/splitted.vcf" > "$outputDir/tmp/$tmpfd/splitted.mtonly.txt"
    	perl -pe 's/^chr//' "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" > "$outputDir/tmp/$tmpfd/leftnormalized2_nochr.vcf"
    	cat "$outputDir/tmp/$tmpfd/leftnormalized2_nochr.vcf" "$outputDir/tmp/$tmpfd/splitted.mtonly.txt" > "$outputDir/tmp/$tmpfd/leftnormalized.vcf"
    else
    grep "^chrMT" "$outputDir/tmp/$tmpfd/splitted.addchr.vcf" > "$outputDir/tmp/$tmpfd/splitted.mtonly.vcf"
    	cat "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" "$outputDir/tmp/$tmpfd/splitted.mtonly.vcf" > "$outputDir/tmp/$tmpfd/leftnormalized.vcf"
    fi
fi
if [ $count == "0" ] && [ $count2 == "0" ]
then
    (bcftools norm -f $genomeRef -o "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" "$outputDir/tmp/$tmpfd/splitted.vcf") 2>&1 | tee -a "$outputDir/tmp/log.txt"
    if [ $count3 == "1" ] 
    then
        perl -pe 's/^([^#])/chr\1/' "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" > "$outputDir/tmp/$tmpfd/leftnormalized.vcf"
    else
    mv "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" "$outputDir/tmp/$tmpfd/leftnormalized.vcf"
    fi
fi;
if [ $count == "1" ] && [ $count2 == "1" ] 
then
    (bcftools norm -f $genomeRef -o "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" "$outputDir/tmp/$tmpfd/splitted.vcf") 2>&1 | tee -a "$outputDir/tmp/log.txt" 
    if [ $count3 == "0" ] 
    then
        perl -pe 's/^chr//' "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" > "$outputDir/tmp/$tmpfd/leftnormalized.vcf"
    else
    mv "$outputDir/tmp/$tmpfd/leftnormalized2.vcf" "$outputDir/tmp/$tmpfd/leftnormalized.vcf"
    fi    
fi 
bgzip -c "$outputDir/tmp/$tmpfd/leftnormalized.vcf" > "$outputDir/tmp/$tmpfd/leftnormalized.vcf.gz"
tabix -p vcf "$outputDir/tmp/$tmpfd/leftnormalized.vcf.gz"
if [[ $cosmicVCF == *.vcf ]]; then bgzip -c "$cosmicVCF" > "$cosmicVCF.gz"; tabix -p vcf "$cosmicVCF.gz"; cosmicVCF="$cosmicVCF.gz"; fi
if [[ $dbSNPVCF == *.vcf ]]; then bgzip -c "$dbSNPVCF" > "$dbSNPVCF.gz"; tabix -p vcf "$dbSNPVCF.gz"; dbSNPVCF="$dbSNPVCF.gz"; fi
## NOTICE: .vcf.gz files will be unzipped, bgzipped and then tabix-ed since it might not be tabix-ed correctly.
if [[ $dbSNPVCF == *.vcf.gz ]]; then tabix -f -p vcf "$dbSNPVCF"; fi
if [[ $cosmicVCF == *.vcf.gz ]]
then
  if tabix -f -p vcf "$cosmicVCF"; then
    echo "The cosmic VCF file can be tabixed."
  else
    cd $(dirname "$cosmicVCF")
    if [ -f "${cosmicVCF::-3}" ]; then rm "${cosmicVCF::-3}"; fi;
    gunzip "$cosmicVCF"
    bgzip -c "${cosmicVCF::-3}" > "$cosmicVCF";
    tabix -f -p vcf "$cosmicVCF"
  fi
fi
### dbsnp annotation via bcftools
bcftools annotate -c ID -a "$dbSNPVCF" "$outputDir/tmp/$tmpfd/leftnormalized.vcf.gz" > "$outputDir/tmp/$tmpfd/dbsnp_anno_tmp.vcf"
count4=$(grep -v '^#' "$outputDir/tmp/$tmpfd/dbsnp_anno_tmp.vcf" | head -1 | grep "^chr" | wc -l)
if [ $count4 == "1" ] 
then
    perl -pe 's/^chr//' "$outputDir/tmp/$tmpfd/dbsnp_anno_tmp.vcf" > "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf"
else
    mv "$outputDir/tmp/$tmpfd/dbsnp_anno_tmp.vcf" "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf"
fi
bgzip -c "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf" > "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf.gz"
tabix -p vcf "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf.gz"
### cosmic annotation via bcftools
bcftools annotate -c ID,+GENE -a "$cosmicVCF" "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf.gz" > "$outputDir/tmp/$tmpfd/cosmic_dbsnp.vcf"
sed '/\trs[0-9]\+\t/d' "$outputDir/tmp/$tmpfd/cosmic_dbsnp.vcf" > "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf"
(java -Xmx4G -jar "$seqtools_snpeff/snpEff.jar" -canon -no-downstream -no-upstream -no-intergenic -no-intron -no-utr -noNextProt -noMotif $genomeVer -s "$outputDir/tmp/$tmpfd/annodbsnpRemove.html" "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf" > "$outputDir/tmp/$tmpfd/snpeff_anno.vcf") 2>&1 | tee -a  "$outputDir/tmp/log.txt"
(cat "$outputDir/tmp/$tmpfd/snpeff_anno.vcf" | java -jar "$seqtools_snpeff/SnpSift.jar" filter "(ANN[*].BIOTYPE = 'protein_coding') | (ANN[*].EFFECT has 'splice')"  > "$outputDir/tmp/$tmpfd/snpeff_proteincoding.vcf") 2>&1 | tee -a  "$outputDir/tmp/log.txt"
# Remove synonymous variants
grep -Ev 'synonymous_variant|start_retained|stop_retained_variant' "$outputDir/tmp/$tmpfd/snpeff_proteincoding.vcf" > "$outputDir/tmp/$tmpfd/nonsyn_splicing.vcf"
grep -wv 'LOW' "$outputDir/tmp/$tmpfd/nonsyn_splicing.vcf" > "$outputDir/tmp/$tmpfd/nonsyn_splicing2.vcf"
read a <<< $(gzip -cd $dbnsfpFile | head -1 | grep 'phyloP7way_vertebrate' | wc -l)
if (( $a > 0 ))
then
dbnsfpField=SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,VEST3_score,CADD_raw,CADD_phred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,GERP++_NR,GERP++_RS,phyloP7way_vertebrate,phastCons7way_vertebrate,SiPhy_29way_logOdds 
else
dbnsfpField=SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,VEST3_score,CADD_raw,CADD_phred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,GERP++_NR,GERP++_RS,phyloP100way_vertebrate,phastCons100way_vertebrate,SiPhy_29way_logOdds 
fi
(java -jar "$seqtools_snpeff/SnpSift.jar" dbNSFP -f $dbnsfpField -v -db "$dbnsfpFile" "$outputDir/tmp/$tmpfd/nonsyn_splicing2.vcf" > "$outputDir/tmp/$tmpfd/nonsyn_splicing_dbnsfp.vcf") 2>&1 | tee -a "$outputDir/tmp/log.txt"
if (( $a > 0 ))
then
(java -jar "$seqtools_snpeff/SnpSift.jar" extractFields -s "," -e "." "$outputDir/tmp/$tmpfd/nonsyn_splicing_dbnsfp.vcf" CHROM POS ID REF ALT "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].GENEID" "ANN[0].FEATURE" "ANN[0].FEATUREID" "ANN[0].BIOTYPE" "ANN[0].HGVS_C" "ANN[0].HGVS_P" dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_LRT_score dbNSFP_LRT_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_FATHMM_score dbNSFP_FATHMM_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_VEST3_score dbNSFP_CADD_raw dbNSFP_CADD_phred dbNSFP_MetaSVM_score dbNSFP_MetaSVM_pred dbNSFP_MetaLR_score dbNSFP_MetaLR_pred dbNSFP_GERP___NR dbNSFP_GERP___RS dbNSFP_phyloP7way_vertebrate dbNSFP_phastCons7way_vertebrate dbNSFP_SiPhy_29way_logOdds > "$outputDir/tmp/$tmpfd/annoTable.txt") 2>&1 | tee -a "$outputDir/tmp/log.txt"
else
(java -jar "$seqtools_snpeff/SnpSift.jar" extractFields -s "," -e "." "$outputDir/tmp/$tmpfd/nonsyn_splicing_dbnsfp.vcf" CHROM POS ID REF ALT "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].GENEID" "ANN[0].FEATURE" "ANN[0].FEATUREID" "ANN[0].BIOTYPE" "ANN[0].HGVS_C" "ANN[0].HGVS_P" dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_LRT_score dbNSFP_LRT_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_FATHMM_score dbNSFP_FATHMM_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_VEST3_score dbNSFP_CADD_raw dbNSFP_CADD_phred dbNSFP_MetaSVM_score dbNSFP_MetaSVM_pred dbNSFP_MetaLR_score dbNSFP_MetaLR_pred dbNSFP_GERP___NR dbNSFP_GERP___RS dbNSFP_phyloP100way_vertebrate dbNSFP_phastCons100way_vertebrate dbNSFP_SiPhy_29way_logOdds > "$outputDir/tmp/$tmpfd/annoTable.txt") 2>&1 | tee -a "$outputDir/tmp/log.txt"
fi
annoTable="$outputFileNameShare"
annoTable+='_annoTable.txt'
sed -r 's/(\[|\])//g' <"$outputDir/tmp/$tmpfd/annoTable.txt" > "$outputDir/tmp/$tmpfd/tmp.txt"
sed -r 's/\#//g' <"$outputDir/tmp/$tmpfd/tmp.txt" > "$outputDir/tmp/$tmpfd/tmp2.txt"
sed -r 's/dbNSFP\_//g' <"$outputDir/tmp/$tmpfd/tmp2.txt" > "$outputDir/tmp/$tmpfd/tmp3.txt"
sed -r 's/ANN0\.//g' <"$outputDir/tmp/$tmpfd/tmp3.txt" > "$outputDir/tmp/$tmpfd/$annoTable"
#rm "$outputDir/tmp/$tmpfd/tmp*.txt"
##run annotation for variants after pre-filtering
(java -Xmx4G -jar "$seqtools_snpeff/snpEff.jar" -canon -noNextProt -noMotif $genomeVer -s "$outputDir/tmp/$tmpfd/annoPrefilter.html" "$outputDir/tmp/$tmpfd/leftnormalized.vcf" > "$outputDir/tmp/$tmpfd/snpEffAnnoAll.vcf") 2>&1 | tee -a  "$outputDir/tmp/log.txt"
(java -jar "$seqtools_snpeff/SnpSift.jar" extractFields -s "," -e "." "$outputDir/tmp/$tmpfd/snpEffAnnoAll.vcf" "ANN[0].EFFECT" > "$outputDir/tmp/$tmpfd/annoTableAll.txt") 2>&1 | tee -a "$outputDir/tmp/log.txt"
chmod +x $RcodeDir/createGeneListSNPEFF.sh
($RcodeDir/./createGeneListSNPEFF.sh "$outputDir/tmp/$tmpfd/$annoTable" "$outputDir/tmp/$tmpfd/annoTableAll.txt" "$outputDir/tmp/$tmpfd") 2>&1 | tee -a "$outputDir/tmp/log.txt"
rm "$outputDir/tmp/$tmpfd/$annoTable"
grep -v "protein_protein_contact" "$outputDir/tmp/$tmpfd/annoTableAfterClean.txt" > "$outputDir/tmp/$tmpfd/$annoTable"
annoVCF="$outputFileNameShare"
annoVCF+='_annotated.vcf'
grep -v "protein_protein_contact" "$outputDir/tmp/$tmpfd/nonsyn_splicing_dbnsfp.vcf" > "$outputDir/tmp/$tmpfd/$annoVCF"
genelistFile="$outputFileNameShare"
genelistFile+='_genelist.txt'
mv "$outputDir/tmp/$tmpfd/genelist.txt" "$outputDir/tmp/$tmpfd/$genelistFile"
read num0 <<< $(grep -v '#' $inputVCF | wc -l)
echo -e ""$inputVCF" contains $num0 lines where each line corresponds to one or multiple variants." >> "$outputDir/tmp/log.txt"
echo -e "There are $num0a variants after pre-filtering." >> "$outputDir/tmp/log.txt"
read num <<< $(grep -v '#' "$outputDir/tmp/$tmpfd/leftnormalized.vcf" | wc -l)
echo -e "There are $num variants after decomposing and normalizing variants.">> "$outputDir/tmp/log.txt"
read num2 <<< $(cut -f 3 "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf" | grep -v '#' | grep 'rs' | wc -l)
echo -e "There are $num2 variants (out of $num variants) that are reported by dbSNP."  >> "$outputDir/tmp/log.txt"
read num3 <<< $(cut -f 3 "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf" | grep -v '#' | grep 'COSM' | wc -l)
echo -e "There are $num3 variants (out of $num variants) that are reported by COSMIC. "  >> "$outputDir/tmp/log.txt"
read num3a <<< $(grep -v '#' "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf" | wc -l)
echo -e "There are $((num3a-num+num2)) variants (out of $num variants) that are reported by both dbSNP and COSMIC VCF files."  >> "$outputDir/tmp/log.txt"
echo -e "After removing these variants reported in dbSNP database and keeping these variants reported in COSMIC database, there are $num3a variants left for analysis."  >> "$outputDir/tmp/log.txt"
read num4 <<< $(grep -v 'CHROM' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l)
read num5 <<< $(grep -v 'Mutation' "$outputDir/tmp/$tmpfd/$genelistFile" | wc -l)
echo -e "There are $num4 variants (out of $num3a variants) that are nonsynonymous or splicing ones, which are kept for further analysis." >> "$outputDir/tmp/log.txt"
read num6 <<< $(grep -v '#' "$outputDir/tmp/$tmpfd/$annoTable" | grep "COSM" | wc -l)
echo -e "There are $num6 variants (out of $num4 variants) that are reported by COSMIC."  >> "$outputDir/tmp/log.txt"
read num_spg <<< $(grep 'splice_region_variant' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_spg variants (out of $num4 variants) that are splicing." >> "$outputDir/tmp/log.txt"
read num_del <<< $(grep 'frameshift_variant' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_del variants (out of $num4 variants) that are frameshift." >> "$outputDir/tmp/log.txt"
read num_sl <<< $(grep 'stop_lost' "$outputDir/tmp/$tmpfd/$annoTable"  | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_sl variants (out of $num4 variants) that are stoploss." >> "$outputDir/tmp/log.txt"
read num_sg <<< $(grep 'stop_gained' "$outputDir/tmp/$tmpfd/$annoTable"  | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_sg variants (out of $num4 variants) that are stopgain." >> "$outputDir/tmp/log.txt"
read num_stl <<< $(grep 'start_lost' "$outputDir/tmp/$tmpfd/$annoTable"  | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_stl variants (out of $num4 variants) that are startloss." >> "$outputDir/tmp/log.txt"
read num_non <<< $(grep 'missense_variant' "$outputDir/tmp/$tmpfd/$annoTable"  | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_non variants (out of $num4 variants) that are mis-sense." >> "$outputDir/tmp/log.txt"
echo -e "There are $num5 genes that are associated with $num4 variants, and the gene list is saved in $outputDir/tmp/$tmpfd/$genelistFile." >> "$outputDir/tmp/log.txt"
### Genenrate a statistics summary table read by python to create a table in a pdf file
if [ -f "$outputDir/tmp/$tmpfd/statistics.txt" ]; then rm "$outputDir/tmp/$tmpfd/statistics.txt"; fi
echo -e "Total number of variants in the raw VCF file $num0" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants left after the filter by QUAL >= $minQual, DP >= $minReadDepth, MQ >= $minMapQual $num0a" >> "$outputDir/tmp/$tmpfd/statistics.txt"
#echo -e "Number of variants after decomposing and left normalization $num" >> "$outputDir/tmp/$tmpfd/statistics.txt"
#echo -e "Number of variants reported in dbSNP database $num2" >> "$outputDir/tmp/$tmpfd/statistics.txt"
#echo -e "Number of variants reported in COSMIC database $num3" >> "$outputDir/tmp/$tmpfd/statistics.txt"
#echo -e "Number of variants reported in both dbSNP and COSMIC database $((num3a-num+num2))" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants remaining after removing variants reported in dbSNP while keeping variants in COSMIC $num3a"  >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants (out of $num3a variants) that are nonsynonymous or splicing ones $num4" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants (out of $num4 variants) that are reported in COSMIC $num6" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of genes associated with $num4 variants $num5" >> "$outputDir/tmp/$tmpfd/statistics.txt"
## Genenrate an effect type summary table read by python to create a table in a pdf file
if [ -f "$outputDir/tmp/$tmpfd/effectType.txt" ]; then rm "$outputDir/tmp/$tmpfd/effectType.txt"; fi
echo -e "Number of variants (out of $num4 variants) that are frameshift $num_del" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are stoploss $num_sl" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are stopgain $num_sg" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are startloss $num_stl" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are mis-sense $num_non" >> "$outputDir/tmp/$tmpfd/effectType.txt"
read num_dins <<< $(grep 'disruptive_inframe_insertion' "$outputDir/tmp/$tmpfd/$annoTable"  | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "Number of variants (out of $num4 variants) that are disruptive inframe insertion $num_dins" >> "$outputDir/tmp/$tmpfd/effectType.txt"
read num_ddel <<< $(grep 'disruptive_inframe_deletion' "$outputDir/tmp/$tmpfd/$annoTable"  | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "Number of variants (out of $num4 variants) that are disruptive inframe deletion $num_ddel" >> "$outputDir/tmp/$tmpfd/effectType.txt"
read num_ins2 <<< $(grep -v 'disruptive_inframe_insertion' "$outputDir/tmp/$tmpfd/$annoTable" | grep 'inframe_insertion' | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "Number of variants (out of $num4 variants) that are inframe insertion $num_ins2" >> "$outputDir/tmp/$tmpfd/effectType.txt"
read num_del2 <<< $(grep -v 'disruptive_inframe_deletion' "$outputDir/tmp/$tmpfd/$annoTable" | grep 'inframe_deletion'  | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "Number of variants (out of $num4 variants) that are inframe deletion $num_del2" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are splicing $num_spg" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Total number of variants that are nonsynonymous or splicing $num4" >> "$outputDir/tmp/$tmpfd/effectType.txt"
## A file saving the directories of outputted files and the raw VCF file
if [ -f "$outputDir/tmp/$tmpfd/fileDetails.txt" ]; then rm "$outputDir/tmp/$tmpfd/fileDetails.txt"; fi
echo -e "$inputVCF" >> "$outputDir/tmp/$tmpfd/fileDetails.txt"
echo -e "$outputDir/$genelistFile" >> "$outputDir/tmp/$tmpfd/fileDetails.txt"
echo -e "$outputDir/$annoTable" >> "$outputDir/tmp/$tmpfd/fileDetails.txt"
echo -e "$outputDir/$annoVCF" >> "$outputDir/tmp/$tmpfd/fileDetails.txt"
## A file saving the filtering parameters
if [ -f "$outputDir/tmp/$tmpfd/filterParam.txt" ]; then rm "$outputDir/tmp/$tmpfd/filterParam.txt"; fi
echo -e "minimum quality $minQual" >> "$outputDir/tmp/$tmpfd/filterParam.txt"
echo -e "minimum read depth $minReadDepth" >> "$outputDir/tmp/$tmpfd/filterParam.txt"
echo -e "minimum mapping quality $minMapQual" >> "$outputDir/tmp/$tmpfd/filterParam.txt"
## run python code to generate summary report
outReportName="$outputFileNameShare"
outReportName+='_summaryReport'
cd "$outputDir/tmp/$tmpfd"
R -e "rmarkdown::render('$RcodeDir/summaryreportSnpEffhtml.Rmd',intermediates_dir='$outputDir/tmp/$tmpfd/',output_dir='$outputDir/tmp/$tmpfd/')" --args "$outputDir/tmp/$tmpfd/"
R -e "rmarkdown::render('$RcodeDir/summaryreportSnpEffpdf.Rmd',intermediates_dir='$outputDir/tmp/$tmpfd/',output_dir='$outputDir/tmp/$tmpfd/')" --args "$outputDir/tmp/$tmpfd/"
cd
if [ $count == "1" ]
then
    perl -pe 's/^([^#])/chr\1/' "$outputDir/tmp/$tmpfd/$annoVCF" > "$outputDir/tmp/$tmpfd/tmp.vcf"
    mv "$outputDir/tmp/$tmpfd/tmp.vcf" "$outputDir/$annoVCF"
fi
if [ $count == "0" ] && [ $count2 == "1" ]
then
    mv "$outputDir/tmp/$tmpfd/$annoVCF" "$outputDir/$annoVCF"
fi
if [ $count == "0" ] && [ $count2 == "0" ]
then
mv "$outputDir/tmp/$tmpfd/$annoVCF" "$outputDir/$annoVCF"
fi
mv "$outputDir/tmp/$tmpfd/summaryreportSnpEffpdf.pdf" "$outputDir/$outReportName.pdf"
mv "$outputDir/tmp/$tmpfd/summaryreportSnpEffhtml.html" "$outputDir/$outReportName.html"
mv "$outputDir/tmp/$tmpfd/$annoTable" "$outputDir/$annoTable"
mv "$outputDir/tmp/$tmpfd/$genelistFile" "$outputDir/$genelistFile"
done
mv "$outputDir/tmp/log.txt" "$outputDir/log.txt"
echo "An annotation table associated with $num4 variants is saved in $outputDir/$annoTable, and a vcf file with annotation information is saved in $outputDir/$annoVCF." >> "$outputDir/log.txt"
rm -r -f "$outputDir/tmp"
if [ -f "$outputDir/inputParameters.txt" ]; then rm "$outputDir/inputParameters.txt"; fi
echo -e "annotationMethod=SnpEff" >> "$outputDir/inputParameters.txt"
echo -e "inputVCFDirectory=$inputVCFDir" >> "$outputDir/inputParameters.txt"
echo -e "inputCosmicVCF=$cosmicVCF" >> "$outputDir/inputParameters.txt"
echo -e "inputDbSNPVCF=$dbSNPVCF" >> "$outputDir/inputParameters.txt"
echo -e "inputGenomeReference=$genomeRef" >> "$outputDir/inputParameters.txt"
echo -e "inputDbNSFPFile=$dbnsfpFile" >> "$outputDir/inputParameters.txt"
echo -e "inputMinVariantQualityScore=$minQual" >> "$outputDir/inputParameters.txt"
echo -e "inputMinReadDepth=$minReadDepth" >> "$outputDir/inputParameters.txt"
echo -e "inputMinMappingQuality=$minMapQual" >> "$outputDir/inputParameters.txt"
echo -e "inputGenomeBuild=$genomeVer" >> "$outputDir/inputParameters.txt"
echo -e "outputDirectory=$outputDir" >> "$outputDir/inputParameters.txt"
if [ -f "$outputDir/errormessage.txt" ]; then rm "$outputDir/errormessage.txt"; fi; 
cnterr1=$(grep -i "fail" "$outputDir/log.txt" | head -1 | wc -l) 
cnterr2=$(grep -i "error" "$outputDir/log.txt" | head -1 | wc -l) 
cnterr3=$(grep -i "Couldn't download the file" "$outputDir/log.txt" | head -1 | wc -l) 
if [ "$cnterr1" != 0 ] || [ "$cnterr2" != 0 ] || [ "$cnterr3" != 0 ] 
then 
grep -i "Couldn't" >> "$outputDir/errormessage.txt" 
grep -i "fail" "$outputDir/log.txt" >> "$outputDir/errormessage.txt" 
grep -i "error" "$outputDir/log.txt" >> "$outputDir/errormessage.txt" 
cd $outputDir 
tar -zcvf "BugReport.tar.gz" "log.txt" "inputParameters.txt" "errormessage.txt" 
cd 
exit 7 
fi

if [ -f "$outputDir/BugReport.tar.gz" ]; then rm "$outputDir/BugReport.tar.gz"; fi; 
if [ ! -f "$outputDir/$outReportName.pdf" ]
then 
cd $outputDir 
tar -zcvf "BugReport.tar.gz" "log.txt" "inputParameters.txt"
cd 
exit 6
fi
echo completed `date +'%Y-%m-%d %T'` >> "$outputDir/log.txt"
