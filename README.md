# MobileAnn
MobileAnn is used to check if SV calls are mobile element insertions. Many SV callers are able to detect mobile element
insertions, but classify these events as other types of variants, such as interchromosomal translocations or inversions.

using a vcf file containing mobile element insertions, MobileAnn searches for these events and changes their classification.
By correcting the classification of these events, interpretation of the variants is made easier.

MobileAnn is  distributed with the AnnotateFrq script, which is used to add frequency and occurance
tags to multip sample vcf files (such as multisample MELT output). Adding these tags, the vcf files may be used as a frequency database. 

# Run

python MobileAnn.py --sv sv.vcf --me me.vcf --rm repeats > updated_sv.vcf

the mobile element vcf (me.vcf) may be produced using MELT or similar callers.

Note: MobileAnn assumes that the naming of the chromosomes are the same within all input files. If you have mixed naming of your input files, you can try to remove the "chr" prefix through the following command:

    sed -i -e 's/chr//g' file

A repeatmaske file may be downloaded via the UCSC tablebrowser.

The AnnotateFrq is  run using the following script

	python AnnotateFrq.py --vcf input.vcf > annotated.vcf

# Install
MobileAnn requires numpy. Numpy may be installed using the following command

	pip install numpy

# Algorithm
    
    MobileAnn accepts three files, a sv vcf file, a mobile element (ME) vcf file, and a repeat tab file. MobileAnn will search the sv file for variants having a breakpoint within a repeat element and a ME call.
    If such variant is found, that SV call is replaced by the ME call.
    If multiple SV calls support a single  ME insertions, the ME insertion will be printed only once.

