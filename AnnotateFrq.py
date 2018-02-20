import argparse

parser = argparse.ArgumentParser("""AnnotateFrq - add SVDB style FRQ and OCC tags to input vcf file""")
parser.add_argument('--vcf', type=str, help="a vcf containing sv", required = True)
parser.add_argument('--frq', type=str, help="frequency tag (default=FRQ)")
parser.add_argument('--occ', type=str, help="occurances tag (default=OCC)")

args = parser.parse_args()


for line in open(args.vcf):
    if line[0] == "#":
        if "CHROM" in line:
            print ("##INFO=<ID={},Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">".format(args.occ))
            print ("##INFO=<ID={},Number=1,Type=Float,Description=\"The frequency of the event in the database\">".format(args.frq))
            print line.strip()
        else:
            print line.strip()
        continue


    
    FRQ=0
    OCC=0
    content=line.strip().split()
    samples=len(content)-9
    for i in range(9,len(content)):
        if not "0/0:" in content[i]:
            OCC+=1
            
    if samples:
        FRQ=OCC/float(samples)
    content[7]+=";{}={};{}={}".format(args.occ,OCC,args.frq,round(FRQ,4))
    print "\t".join(content)
