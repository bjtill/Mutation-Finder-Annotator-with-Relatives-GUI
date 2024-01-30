#!/bin/bash
#B.T. Dec 14, 2023
#GUI version
#V1.1 tests okay
#V1.2 Changes to text descriptions (Jan 24) 
#V1.6 Add fixes directly to 1.2 to avoid confusion. This version tests okay except for no true replicate filtering.
#V1.7 Adds replicate filtering and reporting.  
##################################################################################################################################

zenity --width 1500 --info --title "Mutation Finder and Annotator with Relatives (MFAR) GUI : Click OK to start" --text "
ABOUT: 
This program identifies DNA sequence variants unique to related samples and annotates them using the program SnpEff.  It was built to identify chemically induced point mutations but will work on all biallelic SNVs. The program uses information on biological or technical replicates and genetic relationships (e.g. siblings) to calculate the total number of novel single nucleotide variants (SNVs) in a line.  

This program is a fork of the original Mutation Finder and Annotator (MFA) program that adds multi-sample, and genetic relationship functionality. 
  
PREREQUISITES:
1) A multi-sample VCF filtered for biallelic SNVs that is compressed with bgzip (ending in .vcf.gz), 2) an index of the .vcf.gz file (using tabix), 3) A SnpEff genome database that matches the genome assembly used to generate the VCF file, 4) a three column file listing samples and their relationships. See _ for detailed instructions. 5) the following tools installed: 
SnpEff, SnpSift, Java, bash, zenity, awk, tr, datamash, bcftools. This program was built to run on Ubuntu 20.04 and higher. See the readme file for information on using with other operating systems.  

OUTPUTS:
A summary table is produced that lists the unique mutations in each set of related samples. In addition, SnpEff output files, merged and non-merged VCFs and bcftools stats outputs are produced. Please see the detailed instructions here:  

TO RUN:
Click OK to start. When prompted, enter the name for your analysis directory. A new directory will be created and the files created will be deposited in the directory.  Follow the information to select files and start the program.  


LICENSE:  
MIT License
Copyright (c) 2024 Bradley John Till
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the *Software*), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Version Information:  Version 1.7, January 26, 2024"

directory=`zenity --width 500 --title="DIRECTORY" --text "Enter text to create a new directory (e.g. Sample1234).  
WARNING: No spaces or symbols other than an underscore. DO NOT USE THE SAME NAME AS YOUR SAMPLES FILE!" --entry`

if [ "$?" != 0 ]
then
    exit
    fi
mkdir $directory
cd $directory

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>MFAt.log 2>&1
now=$(date)  
echo "Mutation Finder and Annotator with Relatives (MFAR) GUI, Version 1.7
Script Started $now."  

zenity --width 500 --info --title "VCF Selection" --text "Click OK to select the VCF file. The file should be compressed and end in .vcf.gz "
ZAR=$(zenity --file-selection --title="Select the .vcf.gz file" --file-filter=*.gz)
if [ "$?" != 0 ]
then
    exit
    fi
echo $ZAR >> variantpath

zenity --width 500 --info --title "File Selection" --text "Click OK to select your sample relationship table.  CAUTION: It must be a three column file and samplen names must match the VCF."
samplefile=`zenity --file-selection --title="Select Your sample relationship table."  --save`
cp $samplefile samples.txt

SEG=$(zenity --width 600 --entry --title "Enter the exact name of the SNPeff genome. This is case sensitive." --text="SNPeff Genome Name")
if [ "$?" != 0 ]
then
    exit
    fi
echo $SEG >> SNPeffGenome

zenity --width 500 --info --title "SnpEff Location" --text "Click OK to choose the SnpEff.jar file"
SEP=$(zenity --file-selection --title="Select the SnpEff.jar file")
if [ "$?" != 0 ]
then
    exit
    fi

echo $SEP >> SnpEffpath

zenity --width 500 --info --title "SnpSift Location" --text "Click OK to choose the SnpSift.jar file"
SSP=$(zenity --file-selection --title="Select the SnpSift.jar file")
if [ "$?" != 0 ]
then
    exit
    fi
echo $SSP >> SnpSIFTpath

zenity --width 500 --info --title "READY TO LAUNCH" --text "Click OK to start the Mutation Finder and Annotator with Relatives program. Progress is indicated by a progress bar. A log file titled MFAR.log will be created."

(#Start progress bar
echo "# Finding unique mutations in the selected samples"; sleep 2

##############################################################################################################################################

awk '{print > ($2".reps")}' samples.txt
for i in *.reps; do 
mkdir ${i%.*}
mv $i ./${i%.*}/ 
cp samples.txt ./${i%.*}/ 
cp SnpEffpath ./${i%.*}/
cp variantpath ./${i%.*}/
cp SnpSIFTpath ./${i%.*}/
cp SNPeffGenome ./${i%.*}/
cd ${i%.*}

a=$(wc -l *.reps | awk '{print $1}')
awk -v var=$a 'NR==1 {if (var>1) print "Yes"}' samples.txt > TrueReps

find . -name 'TrueReps' -type f -empty -delete
if [ -f "TrueReps" ]; then

#Assumes GT is first value in saple columns of VCF
a=$(pwd)
b=$(head -1 variantpath)
gunzip -c $b > $a/tmpvcf
grep "#" tmpvcf > header
grep -v "#" tmpvcf > tvcnh
for i in *.reps; do 
awk '{print "awk -v b=\x22"$1"\x22 \x27{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}\x27 tmpvcf >> reps1"}' $i > ${i%.*}.sh
chmod +x ${i%.*}.sh
./${i%.*}.sh
done

awk '{print "awk \x27{print substr($"$1", 1, 3)}\x27 tvcnh > "$1".r2"}' reps1 > reps3.sh 
chmod +x reps3.sh
./reps3.sh #creates GT lists
ls *.r2 > r2list
#Created to work with up to 6 replicates in a relationship set. Note that because the resulting VCF is the input for the snpSift, can remove reference alleles here

a=$(wc -l r2list | awk '{print $1}')
cp r2list ${a}.numrep

awk -v var=$a 'NR==1 {if (var<7) print "Y"}' r2list > ReplicatesStat
find . -name 'ReplicatesStat' -type f -empty -delete
if [ -f "ReplicatesStat" ]; then
###
if [ -f "2.numrep" ]; then
a=$(awk '{print $1".r2"}' reps1 | datamash transpose | tr '\t' ' ' )
paste $a | tr '|' '/' | awk '{if ($1!="0/0" && $1==$2) print NR}' > matchlist
fi

if [ -f "3.numrep" ]; then
a=$(awk '{print $1".r2"}' reps1 | datamash transpose | tr '\t' ' ' )
paste $a | tr '|' '/' | awk '{if ($1!="0/0" && $1==$2 && $2==$3) print NR}' > matchlist
fi

if [ -f "4.numrep" ]; then
a=$(awk '{print $1".r2"}' reps1 | datamash transpose | tr '\t' ' ' )
paste $a | tr '|' '/' | awk '{if ($1!="0/0" && $1==$2 && $2==$3 && $3==$4) print NR}' > matchlist
fi

if [ -f "5.numrep" ]; then
a=$(awk '{print $1".r2"}' reps1 | datamash transpose | tr '\t' ' ' )
paste $a | tr '|' '/' | awk '{if ($1!="0/0" && $1==$2 && $2==$3 && $3==$4 && $4==$5) print NR}' > matchlist
fi

if [ -f "6.numrep" ]; then
a=$(awk '{print $1".r2"}' reps1 | datamash transpose | tr '\t' ' ' )
paste $a | tr '|' '/' | awk '{if ($1!="0/0" && $1==$2 && $2==$3 && $3==$4 && $4==$5 &&) print NR}' > matchlist
fi
Z=`pwd | awk -F/ '{print $NF}'`
awk -v var="$Z" 'NR==1 {print "Acceptable Number of Replicates Selected in Replicate Group", var}' samples.txt > ${Z}.ReplicateWarning
awk 'NR==FNR{data[$1]; next}FNR in data' matchlist tvcnh | cat header - > Reps.vcf

for i in *.reps; do 
awk 'FNR==NR {b[$3]; next}  $3 in b' $i samples.txt > sibs;
awk '{print $1, "1"}' sibs > repcode; done 
grep -vFf sibs samples.txt | awk '{print $1, "2"}' > nonrelated
cat repcode nonrelated > SSlist
awk '{print > ($1".b")}' SSlist
for i in *.b; do 
awk '{if ($2==2) print "isRef (GEN["$1"]) &"; else if ($2==1) print "(isHom (GEN["$1"]) & isVariant(GEN["$1"]) || isHet ( GEN["$1"])) & "}' $i > ${i%.*}.bb
done

paste *.bb > script.stage
#Clean up the tabs and get rid of terminal &
tr '\t' ' ' < script.stage > script.tab
rev script.tab | cut -c3- | rev > script2.tab
#Collect the paths for SnpSift and the VCF file and create a shell script to run SnpSift
a=$(head -1 SnpSIFTpath)
awk -v awkvar1="$a" '{print "#!/bin/bash" "\n" "#Induced mutation identification" "\n" "java -Xmx64g -jar " awkvar1 " filter \x22"$0"\x22 Reps.vcf > output.vcf"}' script2.tab > MutHunt.sh

# give permission for the shell script to run.  
chmod +x MutHunt.sh

else 
Z=`pwd | awk -F/ '{print $NF}'`
awk -v var="$Z" 'NR==1 {print "ERROR: More than 6 related replicates selected in replicate group", var}' samples.txt > ${Z}.ReplicateWarning

fi 

else
#############################################End handling true reps##########################################################



Z=`pwd | awk -F/ '{print $NF}'`
awk -v var="$Z" 'NR==1 {print "No replicates identified in replicate group", var}' samples.txt > ${Z}.ReplicateWarning

for i in *.reps; do 
awk 'FNR==NR {b[$3]; next}  $3 in b' $i samples.txt > sibs;
awk '{print $1, "1"}' $i > repcode; done 
grep -vFf sibs samples.txt | awk '{print $1, "2"}' > nonrelated
cat repcode nonrelated > SSlist
awk '{print > ($1".b")}' SSlist
for i in *.b; do 
awk '{if ($2==2) print "isRef (GEN["$1"]) &"; else if ($2==1) print "(isHom (GEN["$1"]) & isVariant(GEN["$1"]) || isHet ( GEN["$1"])) & "}' $i > ${i%.*}.bb
done

paste *.bb > script.stage
#Clean up the tabs and get rid of terminal &
tr '\t' ' ' < script.stage > script.tab
rev script.tab | cut -c3- | rev > script2.tab
#Collect the paths for SnpSift and the VCF file and create a shell script to run SnpSift
a=$(head -1 SnpSIFTpath)
b=$(head -1 variantpath)
awk -v awkvar1="$a" -v awkvar2="$b" '{print "#!/bin/bash" "\n" "#Induced mutation identification" "\n" "java -Xmx64g -jar " awkvar1 " filter \x22"$0"\x22 " awkvar2 " > output.vcf"}' script2.tab > MutHunt.sh

# give permission for the shell script to run.  
chmod +x MutHunt.sh
fi
echo "75" 
echo "# Running genotype selection script"; sleep 2 
cp *ReplicateWarning ..
./MutHunt.sh

echo "85" 
echo "# Running SNPeff on resulting VCF"; sleep 2
 
#Collect the path for SnpEff and the name of the genome database to be used for SnpEff

#Run SnpEff				

c=$(head -1 SnpEffpath)
d=$(head -1 SNPeffGenome)

java -Xmx32g -jar $c $d output.vcf > annotated_ouput.vcf

a=$(head -1 SSlist | awk '{print $1}')
mv snpEff_genes.txt ${a}_snpEff_genes.txt
mv snpEff_summary.html ${a}_snpEff_summary.html
mv ${a}_snpEff_genes.txt ..
mv ${a}_snpEff_summary.html ..

awk '{print "cp annotated_ouput.vcf", $1"_annotated_output.vcf"}' *.reps > copy.sh
chmod +x copy.sh 
./copy.sh 
cp *annotated_output.vcf ..
rm *.b *.bb *.sh nonrelated repcode samples.txt script.stage script.tab script2.tab sibs SNPeffGenome SnpEffpath SnpSIFTpath SSlist variantpath 
echo "90" 
cd .. 
done 
for i in *.vcf; do 
bgzip $i 
bcftools index ${i%.*}.vcf.gz
done 
echo "# Counting Variants and Cleaning Up"; sleep 2
awk '{print > ($3".sibs")}' samples.txt
for i in *.sibs; do 
awk '{print $1"_annotated_output.vcf.gz"}' $i > ${i%.*}.merge; done 
for i in *.merge; do 
a=$(sed 's/_annotated_output.vcf.gz//g' $i | datamash transpose | tr '\t' '_')
bcftools concat -f $i -D -a >  ${a}_merge.vcf; done 
for i in *_merge.vcf; do 
bcftools stats $i > ${i%_*}.stats; done 
#Collect stats for each line.

for i in *.stats; do  awk '{if ($5 == "records:") print FILENAME, $6}' $i | sed 's/.stats//g' > ${i%.*}.st2; done 

now=$(date +"%m_%d_%Y")
cat *.st2 | awk 'BEGIN{print "Related_Samples", "Total_Unique_Biallelic_SNVs"}1'> RelativesVariantCount_${now}.text


echo "99" 


echo "# Final processing steps.  Program almost finished."; sleep 2
) | zenity --width 800 --title "PROGRESS" --progress --auto-close
now=$(date)  
echo "Program finished $now."
#Collect information to add to log file.  The program is not technically finished, but the above line ensures an end time appears on the log file. 

awk '{print "SNPeff genome used:", $1}' SNPeffGenome > sng
awk '{print "Path to VCF file used:", $1}' variantpath > vpath
awk '{print "Path to SnpEFF.jar:", $1}' SnpEffpath > sepath
awk '{print "Path to SnpSIFT.jar:", $1}' SnpSIFTpath > sspath
echo "Samples selected with replicate and relatives relationships:" > sibrelation
echo " " > space
cat *ReplicateWarning > Repstat
echo "Status of Replicates:" > Rephead
cat MFAt.log space sng space vpath space sspath space sepath space sibrelation samples.txt space Rephead Repstat > MFAR.log
#Remove temporary files 
rm *.b *.bb MutHunt.sh script.stage script.tab script2.tab samplelist1 SL5 SSlist mut1 mut2 Mutants MFAt.log variantpath vpath SnpEffpath sepath SnpSIFTpath sspath sng 
rm MutantL MutantTab SNPeffGenome sibrelation *.merge *.st2 space
#get rid of empty vcf files
find . -name '*.vcf' -type f -empty -delete


#Cleaning
mkdir ForTrash
awk '{print "mv", $2, "./ForTrash/"}' samples.txt > mover.sh
chmod +x mover.sh 
./mover.sh
mkdir BCFtools_Stats
mv *.stats ./BCFtools_Stats/
mkdir Merged_AnnotatedVCFs
mv *.vcf ./Merged_AnnotatedVCFs/
rm *.sibs
mkdir NonMerged_AnnotatedVCFs  #If you don't have either relatives or tech reps
mv *.gz ./NonMerged_AnnotatedVCFs/
mv *.csi ./NonMerged_AnnotatedVCFs/
rm samples.txt mover.sh *.ReplicateWarning Rephead Repstat
mkdir snpEff_genes_summary
mv *.html ./snpEff_genes_summary/
mv *_genes.txt ./snpEff_genes_summary/
#NOTE: snpEff.html and _genes.txt are for the single sample prior to relatives merging

#End of program 
##################################################################################################################################


