# Mutation-Finder-Annotator-with-Relatives-GUI
A tool with graphical user interface to facilitate finding and annotating induced mutations. This is a modification of the MFA tool designed to work on multiple samples and genetically related samples. 

______________________________________________________________________________________________________________________________________________

Use at your own risk.
I cannot provide support. All information obtained/inferred with this script is without any implied warranty of fitness for any purpose or use whatsoever.

ABOUT: 
This program identifies DNA sequence variants unique to related samples and annotates them using the program SnpEff.  It was built to identify chemically induced point mutations but will work on all biallelic SNVs. The program uses information on biological or technical replicates and genetic relationships (e.g. siblings) to calculate the total number of novel single nucleotide variants (SNVs) in a line.  

This program is a fork of the original Mutation Finder and Annotator (MFA) program that adds multi-sample, and genetic relationship functionality.  See https://github.com/bjtill/Mutation-Finder-Annotator-CLI, and associated links and publication for more details.  Detailed instructions are included in the MFAR_Example_Data directory.  
  
PREREQUISITES:
1) A multi-sample VCF filtered for biallelic SNVs that is compressed with bgzip (ending in .vcf.gz).
2) An index of the .vcf.gz file (using tabix).
3) A SnpEff genome database that matches the genome assembly used to generate the VCF file.
4) A three column file listing samples and their relationships.  Column 1 contains the sample name exactly as it appears in the VCF file.  Column 2 contains the replicate status.  Capital letters are used starting with A.  All samples are assigned a different letter except for technical or biological replicates which are assigned the same letter. The third column designates relationship.  Samples that derive from the same mutant line are assigned the same letter (see MFAR_Example_Data directory for instructions).
 5) The following tools installed: SnpEff, SnpSift, Java, bash, zenity, awk, tr, datamash, bcftools. This program was built to run on Ubuntu 20.04 and higher. See https://github.com/bjtill/Mutation-Finder-Annotator-GUI for information on using with other operating systems.  

OUTPUTS:
A summary table is produced that lists the unique mutations in each set of related samples. In addition, SnpEff output files, merged and non-merged VCFs and bcftools stats outputs are produced. Please see the detailed instructions in the  MFAR_Example_Data directory. 

TO RUN:
Click OK to start. When prompted, enter the name for your analysis directory. A new directory will be created and the files created will be deposited in the directory.  Follow the information to select files and start the program.  
