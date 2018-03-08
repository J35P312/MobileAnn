# MobileAnn
MobileAnn is used to add frequency annotation to Mobile Element (ME) Vcf files, and to reclassify SV vcf calls into mobile element calls.
# Run
MobileAnn consists of three modules. The SV annotation module, ME annotatation module, and the frequency annotation module.

The SV annotatation module compares SV calls to mobile element calls, and changes the SV calls into ME calls. The purpose of this module is to re-classify SV calls that represent ME calls: Many SV callers are able to detect Mobile elements, but will classify these events as translocations or inversions. The SV-annotation module is run through the following command:


    python MobileAnn.py --sv_annotate --sv_ --sv sv.vcf --me me.vcf --rm repeats > updated_sv.vcf

the mobile element vcf (me.vcf) may be produced using MELT or similar callers.

Note: MobileAnn assumes that the naming of the chromosomes are the same within all input files. If you have mixed naming of your input files, you can try to remove the "chr" prefix through the following command:

    sed -i -e 's/chr//g' file

A repeatmasker file may be downloaded via the UCSC tablebrowser.

The frequency annotation module is run using the following command:

	python MobileAnn.py --frq --vcf input.vcf > annotated.vcf

This  module will generate a new vcf file containing the frequency and counts of each mobile element found in the input vcf. The input vcf is asssumed to be a multi-individual vcf, and  the frequencies wiwll be based on the individuals on that vcf.

Lastly the ME_annotation module is run through the following command:

    python MobileAnn.py --me_annotate --vcf <me.vcf> --medb <multi_sample_ME.vcf> > annotated.vcf

THe me_db file is a multi-sample ME file produced through the frequency annotation module. And  the vcf file is a vcf file produced by MELT. THE vcf file is annotated based on the frequencices of the me_db file. 
A mobile element in the medb file, and a mobile element in the vcf file will be considered the same if they have the same alt-entry (such as INS:ALU), and if they are located within the --radius distance (defualt = 50 bp)

# Install
MobileAnn requires numpy. Numpy may be installed using the following command

	pip install numpy

# Algorithm
    
    MobileAnn accepts three files, a sv vcf file, a mobile element (ME) vcf file, and a repeat tab file. MobileAnn will search the sv file for variants having a breakpoint within a repeat element and a ME call.
    If such variant is found, that SV call is replaced by the ME call.
    If multiple SV calls support a single  ME insertions, the ME insertion will be printed only once.

