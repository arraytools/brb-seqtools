#!/bin/bash
seqtools_annovar=/home/brb/annovar
seqtools_annovardata=/home/brb/annovar/humandb
export bdge_bcftools_PATH=/opt/SeqTools/bin/bcftools-1.3
export PATH=$bdge_bcftools_PATH:$PATH

export bdge_htslib_PATH=/opt/SeqTools/bin/samtools-1.3/htslib-1.3
export PATH=$bdge_htslib_PATH:$PATH

inputVCFDir="/home/brb/testdata/GSE48215subset/output"
cosmicVCF="/home/brb/variantAnnoDatabase/cosmic/GRCh37/CosmicCodingMuts.vcf.gz"
dbSNPVCF="/home/brb/variantAnnoDatabase/dbsnp/GRCh37/common_all_20160601.vcf.gz"
genomeRef="/home/brb/testdata/hg19/chr1.fa"
minQual=20
minReadDepth=5
minMapQual=1
genomeVer=hg19
chmod +x $seqtools_annovar/*.pl
RcodeDir="$(dirname /home/brb/seqtools-dl/SeqTools)/code"
if [ ! -d "$RcodeDir" ]; then exit 8; fi;outputDir="/home/brb/testdata/GSE48215subset/output/annovar"
if [ ! -d "$outputDir" ]; then mkdir "$outputDir"; fi
if ls "$outputDir"/*.txt >/dev/null; then rm "$outputDir"/*.txt ; fi
if ls "$outputDir"/*_annotated.vcf >/dev/null; then rm "$outputDir"/*_annotated.vcf ; fi
if ls "$outputDir"/*.pdf >/dev/null; then rm "$outputDir"/*.pdf ; fi
if ls "$outputDir"/*.html >/dev/null; then rm "$outputDir"/*.html ; fi
if ls "$outputDir"/run_variantannose.sh >/dev/null; then rm "$outputDir"/run_variantannose.sh ; fi
if [ -d "$outputDir/tmp" ]; then rm -r -f "$outputDir/tmp";  fi
mkdir "$outputDir/tmp"
if [ -f "$outputDir/tmp/log.txt" ]; then rm "$outputDir/tmp/log.txt"; fi
echo started `date +'%Y-%m-%d %T'` >> "$outputDir/tmp/log.txt"
count2=$(head -1 "$genomeRef" | grep "chr" | wc -l)
if [[ $dbSNPVCF == *.vcf ]]; then count3=$(grep -v '^#' "$dbSNPVCF" | head -1 | grep "chr" | wc -l); fi;
if [[ $dbSNPVCF == *.vcf.gz ]]; then 
    count3=$(zgrep "^chr[0-9]" "$dbSNPVCF" | head -1 | wc -l);
fi;
for inputVCF in $inputVCFDir/*.vcf
do
tmpname="$(basename $inputVCF)"
outputFileNameShare="${tmpname::-4}"
tmpfd="tmp_$outputFileNameShare"
if [ ! -d "$outputDir/tmp/$tmpfd" ]; then mkdir "$outputDir/tmp/$tmpfd";  fi

echo "running sample: $inputVCF" >> "$outputDir/tmp/log.txt"
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
#bgzip -c "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf" > "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf.gz"
#tabix -p vcf "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf.gz"
##download database files required for running ANNOVAR annotation
file1=$genomeVer
file1+="_refGene.txt"
file2=$genomeVer
file2+="_refGeneMrna.fa"
file3=$genomeVer
file3+="_refGeneVersion.txt"
if [ ! -f "$seqtools_annovardata/$file1" ] || [ ! -f "$seqtools_annovardata/$file2" ] || [ ! -f "$seqtools_annovardata/$file3" ] ; then
    (perl $seqtools_annovar/annotate_variation.pl -buildver $genomeVer -downdb -webfrom annovar refGene "$seqtools_annovardata/") 2>&1 | tee -a  "$outputDir/tmp/log.txt";
fi;
file1=$genomeVer
file1+="_knownGene.txt"
file2=$genomeVer
file2+="_knownGeneMrna.fa"
if [ ! -f "$seqtools_annovardata/$file1" ] || [ ! -f "$seqtools_annovardata/$file2" ] ; then
    (perl $seqtools_annovar/annotate_variation.pl -buildver $genomeVer -downdb -webfrom annovar knownGene "$seqtools_annovardata/") 2>&1 | tee -a  "$outputDir/tmp/log.txt";
fi
if [ "$genomeVer" != "hg38" ] ; then
  file1=$genomeVer
  file1+="_ensGene.txt"
  file2=$genomeVer
  file2+="_ensGeneMrna.fa"
  if [ ! -f "$seqtools_annovardata/$file1" ] || [ ! -f "$seqtools_annovardata/$file2" ] ; then
    (perl $seqtools_annovar/annotate_variation.pl -buildver $genomeVer -downdb -webfrom annovar ensGene "$seqtools_annovardata/") 2>&1 | tee -a  "$outputDir/tmp/log.txt";
  fi
fi
file1=$genomeVer
file1+="_dbnsfp30a.txt"
file2=$genomeVer
file2+="_dbnsfp30a.txt.idx"
if [ ! -f "$seqtools_annovardata/$file1" ] || [ ! -f "$seqtools_annovardata/$file2" ] ; then
    (perl $seqtools_annovar/annotate_variation.pl -buildver $genomeVer -downdb -webfrom annovar dbnsfp30a "$seqtools_annovardata/") 2>&1 | tee -a  "$outputDir/tmp/log.txt";
fi
## convert vcf file to avinput, a format used by annovar
(perl $seqtools_annovar/convert2annovar.pl -format vcf4 "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf" > "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.avinput" -includeinfo -comment) 2>&1 | tee -a  "$outputDir/tmp/log.txt"
## identify splicing variants and nonsynonymous SNPs
(perl $seqtools_annovar/variants_reduction.pl "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.avinput" "$seqtools_annovardata/" -protocol nonsyn_splicing,dominant -operation g,m -out "$outputDir/tmp/$tmpfd/reduced" -buildver $genomeVer) 2>&1 | tee -a  "$outputDir/tmp/log.txt"
## generate a vcf file from the reduced results in .avinput
cut -f8- "$outputDir/tmp/$tmpfd/reduced.step2.varlist" > "$outputDir/tmp/$tmpfd/tmpfile"
grep '#' "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf" > "$outputDir/tmp/$tmpfd/metainfo"
cat "$outputDir/tmp/$tmpfd/metainfo" "$outputDir/tmp/$tmpfd/tmpfile" > "$outputDir/tmp/$tmpfd/reduced.vcf"
rm "/$outputDir/tmp/$tmpfd/tmpfile"
rm "/$outputDir/tmp/$tmpfd/metainfo"
## annotate the reduced vcf through annovar
perl "$seqtools_annovar/convert2annovar.pl" -includeinfo -allsample -withfreq -format vcf4 "$outputDir/tmp/$tmpfd/reduced.vcf" > "$outputDir/tmp/$tmpfd/test.avinput"
(perl "$seqtools_annovar/table_annovar.pl" "$outputDir/tmp/$tmpfd/reduced.vcf" "$seqtools_annovardata/" -buildver $genomeVer -out "$outputDir/tmp/$tmpfd/$outputFileNameShare" -remove -protocol refGene,dbnsfp30a -operation g,f -nastring . -vcfinput) 2>&1 | tee -a "$outputDir/tmp/log.txt"
annoTable="$outputFileNameShare"
annoTable+='_annoTable.txt'
annoVCF="$outputFileNameShare"
annoVCF+='_annotated.vcf'
file="$outputFileNameShare.$genomeVer"
file+='_multianno'
mv "$outputDir/tmp/$tmpfd/$file.txt" "$outputDir/tmp/$tmpfd/$annoTable"
mv "$outputDir/tmp/$tmpfd/$file.vcf" "$outputDir/tmp/$tmpfd/$annoVCF"
## re-arrange the annotation table
cut -f -7,9-44,50 "$outputDir/tmp/$tmpfd/$annoTable" > "$outputDir/tmp/$tmpfd/tmp0.txt"
grep -v 'Chr' "$outputDir/tmp/$tmpfd/tmp0.txt" > "$outputDir/tmp/$tmpfd/tmp.txt"
grep 'Chr' "$outputDir/tmp/$tmpfd/tmp0.txt" > "$outputDir/tmp/$tmpfd/tmph.txt"
awk '{ printf("%s\tCOSMIC ID\n", $0); }' "$outputDir/tmp/$tmpfd/tmph.txt" > "$outputDir/tmp/$tmpfd/tmph1.txt"
cat "$outputDir/tmp/$tmpfd/tmph1.txt" "$outputDir/tmp/$tmpfd/tmp.txt" > "$outputDir/tmp/$tmpfd/tmp0.txt"
cut -f 1-5 "$outputDir/tmp/$tmpfd/tmp0.txt" > "$outputDir/tmp/$tmpfd/tmp1.txt"
cut -f 6-43 "$outputDir/tmp/$tmpfd/tmp0.txt" > "$outputDir/tmp/$tmpfd/tmp2.txt"
cut -f 44 "$outputDir/tmp/$tmpfd/tmp0.txt" > "$outputDir/tmp/$tmpfd/tmp3.txt"
paste -d'\t' "$outputDir/tmp/$tmpfd/tmp1.txt" <(cat "$outputDir/tmp/$tmpfd/tmp3.txt") > "$outputDir/tmp/$tmpfd/tmp4.txt"
paste -d'\t' "$outputDir/tmp/$tmpfd/tmp4.txt" <(cat "$outputDir/tmp/$tmpfd/tmp2.txt") > "$outputDir/tmp/$tmpfd/$annoTable"
rm "$outputDir/tmp/$tmpfd/tmp*.txt"
## reorder the genelist file and save it as genelist.txt
grep 'Number_of_deleterious_alleles' "$outputDir/tmp/$tmpfd/reduced.step2.genelist" > "$outputDir/tmp/$tmpfd/tmpa"
grep -v 'Number_of_deleterious_alleles' "$outputDir/tmp/$tmpfd/reduced.step2.genelist"  | sort > "$outputDir/tmp/$tmpfd/tmpb"
genelistFile="$outputFileNameShare"
genelistFile+='_genelist.txt'
cat "$outputDir/tmp/$tmpfd/tmpa" "$outputDir/tmp/$tmpfd/tmpb" > "$outputDir/tmp/$tmpfd/$genelistFile"
cut -f1,3 "$outputDir/tmp/$tmpfd/$genelistFile" > "$outputDir/tmp/$tmpfd/tmp_gl.txt"

rm "/$outputDir/tmp/$tmpfd/tmpa"
rm "/$outputDir/tmp/$tmpfd/tmpb"
## gene annotation through refseq, ensembl and ucsc, for the purpose of drawing figures in R; for hg38, ensembl annotation was not provided
if [ "$genomeVer" == "hg38" ]
  then
	perl $seqtools_annovar/table_annovar.pl "$outputDir/tmp/$tmpfd/leftnormalized.vcf" "$seqtools_annovardata/" -buildver $genomeVer -out "$outputDir/tmp/$tmpfd/anno_all" -remove -protocol refGene,knownGene -operation g,g -nastring . -vcfinput >> "$outputDir/tmp/log.txt"
	file=anno_all.$genomeVer
  else
	perl $seqtools_annovar/table_annovar.pl "$outputDir/tmp/$tmpfd/leftnormalized.vcf" "$seqtools_annovardata/" -buildver $genomeVer -out "$outputDir/tmp/$tmpfd/anno_all" -remove -protocol refGene,knownGene,ensGene -operation g,g,g -nastring . -vcfinput >> "$outputDir/tmp/log.txt"
	file=anno_all.$genomeVer
fi	
file+='_multianno'
if [ "$genomeVer" == "hg38" ]
  then
	cut -d$'\t' -f1-15 "$outputDir/tmp/$tmpfd/$file.txt" > "$outputDir/tmp/$tmpfd/anno_all.txt"
  else
	cut -d$'\t' -f1-20 "$outputDir/tmp/$tmpfd/$file.txt" > "$outputDir/tmp/$tmpfd/anno_all.txt"
fi
## run R code to generate figures
chmod +x $RcodeDir/statsVariantGeneAnno_Auto.sh
($RcodeDir/./statsVariantGeneAnno_Auto.sh "$outputDir/tmp/$tmpfd/anno_all.txt" "$outputDir/tmp/$tmpfd") 2>&1 | tee -a "$outputDir/tmp/log.txt"
### Generate a log file
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
echo -e "There are $((num3a-num+num2)) variants (out of $num variants) that are reported by both dbSNP build 144 and COSMIC v74."  >> "$outputDir/tmp/log.txt"
echo -e "After removing these variants reported in dbSNP database and keeping these variants reported in COSMIC database, there are $num3a variants left for analysis."  >> "$outputDir/tmp/log.txt"
read num4 <<< $(wc -l "$outputDir/tmp/$tmpfd/reduced.step2.varlist" | awk '{print $1}')
read num5 <<< $(grep -v 'Number_of_deleterious_alleles' "$outputDir/tmp/$tmpfd/$genelistFile" | wc -l)
echo -e "There are $num4 variants (out of $num3a variants) that are nonsynonymous or splicing ones, which are kept for further analysis." >> "$outputDir/tmp/log.txt"
read num6 <<< $(grep -v '#' "$outputDir/tmp/$tmpfd/reduced.vcf" | grep "COSM" | wc -l)
echo -e "There are $num6 variants (out of $num4 variants) that are reported by COSMIC v74."  >> "$outputDir/tmp/log.txt"
read num_spg <<< $(grep -v 'exonic;splicing' "$outputDir/tmp/$tmpfd/$annoTable" | grep 'splicing' | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_spg variants (out of $num4 variants) that are splicing." >> "$outputDir/tmp/log.txt"
read num_del <<< $(grep 'frameshift deletion' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_del variants (out of $num4 variants) that are frameshift deletion." >> "$outputDir/tmp/log.txt"
read num_ins <<< $(grep 'frameshift insertion' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_ins variants (out of $num4 variants) that are frameshift insersion." >> "$outputDir/tmp/log.txt"
read num_sl <<< $(grep 'stoploss' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_sl variants (out of $num4 variants) that are stoploss." >> "$outputDir/tmp/log.txt"
read num_sg <<< $(grep 'stopgain' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_sg variants (out of $num4 variants) that are stopgain." >> "$outputDir/tmp/log.txt"
read num_non <<< $(grep 'nonsynonymous' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_non variants (out of $num4 variants) that are nonsynonymous." >> "$outputDir/tmp/log.txt"
read num_spg <<< $((num4-num_del-num_ins-num_sl-num_sg-num_non)) >> "$outputDir/tmp/log.txt"echo -e "There are $num_spg variants (out of $num4 variants) that are splicing." >> "$outputDir/tmp/log.txt"echo -e "There are $num5 genes that are associated with $num4 variants, and the gene list is saved in "$outputDir/tmp/$tmpfd/$genelistFile"." >> "$outputDir/tmp/log.txt"
### Genenrate a statistics summary table read by python to create a table in a pdf file
if [ -f "$outputDir/tmp/$tmpfd/statistics.txt" ]; then rm "$outputDir/tmp/$tmpfd/statistics.txt"; fi
echo -e "Total number of variants in the raw VCF file $num0" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants left after the filter QUAL >= $minQual, DP >= $minReadDepth, MQ >= $minMapQual $num0a" >> "$outputDir/tmp/$tmpfd/statistics.txt"
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
echo -e "Number of variants (out of $num4 variants) that are frameshift deletion $num_del" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are frameshift insertion $num_ins" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are stoploss $num_sl" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are stopgain $num_sg" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are mis-sense $num_non" >> "$outputDir/tmp/$tmpfd/effectType.txt"
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
## A file saving the genome build info
if [ -f "$outputDir/tmp/$tmpfd/genomeBuild.txt" ]; then rm "$outputDir/tmp/$tmpfd/genomeBuild.txt"; fi
echo -e "$genomeVer" >> "$outputDir/tmp/$tmpfd/genomeBuild.txt"
## run python code to generate summary report
outReportName="$outputFileNameShare"
outReportName+='_summaryReport'
cd "$outputDir/tmp/$tmpfd"
R -e "rmarkdown::render('$RcodeDir/summaryreportANNOVARhtml.Rmd',intermediates_dir='$outputDir/tmp/$tmpfd/',output_dir='$outputDir/tmp/$tmpfd/')" --args "$outputDir/tmp/$tmpfd/"
R -e "rmarkdown::render('$RcodeDir/summaryreportANNOVARpdf.Rmd',intermediates_dir='$outputDir/tmp/$tmpfd/',output_dir='$outputDir/tmp/$tmpfd/')" --args "$outputDir/tmp/$tmpfd/"
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
mv "$outputDir/tmp/$tmpfd/summaryreportANNOVARpdf.pdf" "$outputDir/$outReportName.pdf"
mv "$outputDir/tmp/$tmpfd/summaryreportANNOVARhtml.html" "$outputDir/$outReportName.html"
mv "$outputDir/tmp/$tmpfd/$annoTable" "$outputDir/$annoTable"
mv "$outputDir/tmp/$tmpfd/tmp_gl.txt" "$outputDir/$genelistFile"
done
mv "$outputDir/tmp/log.txt" "$outputDir/log.txt"
echo "An annotation table associated with $num4 variants is saved in $outputDir/$annoTable, and a vcf file with annotation information is saved in $outputDir/$annoVCF." >> "$outputDir/log.txt"
rm -r -f "$outputDir/tmp"
if [ -f "$outputDir/inputParameters.txt" ]; then rm "$outputDir/inputParameters.txt"; fi
echo -e "annotationMethod=ANNOVAR" >> "$outputDir/inputParameters.txt"
echo -e "inputVCFDirectory=$inputVCFDir" >> "$outputDir/inputParameters.txt"
echo -e "inputCosmicVCF=$cosmicVCF" >> "$outputDir/inputParameters.txt"
echo -e "inputDbSNPVCF=$dbSNPVCF" >> "$outputDir/inputParameters.txt"
echo -e "inputGenomeReference=$genomeRef" >> "$outputDir/inputParameters.txt"
echo -e "inputMinVariantQualityScore=$minQual" >> "$outputDir/inputParameters.txt"
echo -e "inputMinReadDepth=$minReadDepth" >> "$outputDir/inputParameters.txt"
echo -e "inputMinMappingQuality=$minMapQual" >> "$outputDir/inputParameters.txt"
echo -e "inputGenomeBuild=$genomeVer" >> "$outputDir/inputParameters.txt"
echo -e "outputDirectory=$outputDir" >> "$outputDir/inputParameters.txt"
if [ -f "$outputDir/errormessage.txt" ]; then rm "$outputDir/errormessage.txt"; fi; 
cnterr1=$(grep -i "fail" "$outputDir/log.txt" | head -1 | wc -l) 
cnterr2=$(grep -i "error" "$outputDir/log.txt" | head -1 | wc -l) 
if [ "$cnterr1" != 0 ] || [ "$cnterr2" != 0 ] 
then 
grep -i "fail" "$outputDir/log.txt" >> "$outputDir/errormessage.txt" 
grep -i "error" "$outputDir/log.txt" >> "$outputDir/errormessage.txt" 
cd $outputDir 
tar -zcvf "BugReport.tar.gz" "log.txt" "inputParameters.txt" "errormessage.txt" 
cd 
exit 7 
fi

if [ -f "$outputDir/BugReport.tar.gz" ]; then rm "$outputDir/error.tar.gz"; fi; 
if [ ! -f "$outputDir/$outReportName.pdf" ]
then 
cd $outputDir 
tar -zcvf "BugReport.tar.gz" "log.txt" "inputParameters.txt"
cd 
fi
echo completed `date +'%Y-%m-%d %T'` >> "$outputDir/log.txt"
