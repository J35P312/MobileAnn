import sys
import numpy
import readVCF
import argparse
import sqlite3
import time

def merge_info(SV_line,ME_line):
    merged={}
    for entry in ME_line.split(";"):
        merged[entry.split("=")[0]]=entry.split("=")[1]

    for entry in SV_line.split(";"):
        entry_id=entry.split("=")[0]
        entry_content=entry.split("=")[1]

        if not entry_id in merged and not entry_id in ["END","CIPOS","CIEND","MATEID","SVTYPE","SVLEN"]:
            merged[entry_id]=entry_content

    merged_info=[]
    for entry in merged:
            merged_info.append("{}={}".format(entry,merged[entry]))

    return( ";".join(merged_info) )

def construct_header(args):
    header={};
    header["ALT"]={}
    header["INFO"]={}
    header["FILTER"]={}
    header["FORMAT"]={}
    header["CONTIGS"]=[]
    contigs=False
    reference=""
    subheader={}
    columns=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
    first=True
    first_vcf_header=""
    columns=""
    for vcf in [args.sv,args.me]:
        for line in open(vcf):
            if first:
                if "source=" in line:
                    print line.strip()
                    continue
                if "##file" in line:
                    print line.strip()
                    continue
                if("#CHROM\tPOS" in line):
                    columns=line.strip()

            if(line[0] == "#"):

                if line[0] == line[1] and "=" in line:
                    if("ID=" in line and not "##contig=<ID=" in line):
                        field=line.split("=")[2].split(",")[0]
                        key= line.strip("#").split("=")[0]
                        if not key in header:
                            header[key]={}
                        header[key][field]=line
                    elif "##contig=<ID=" in line and not contigs:
                        header["CONTIGS"].append(line)
                    elif "##reference=" in line and not contigs:
                        reference=line
                    elif not "##source=" in line and not "##file" in line and not "##reference=" in line and not "##contig" in line:
                        key=line.strip("#").split("=")[0]
                        if not key in subheader:
                                subheader[key]=line
            else:
                if header["CONTIGS"]:
                    contigs = True
                break

        first=False

    #print the mandatory header lines in the correct order
    for entry in sorted(header["ALT"]):
        print(header["ALT"][entry].strip())
    del header["ALT"]
    for entry in sorted(header["INFO"]):
        print(header["INFO"][entry].strip())
    del header["INFO"]
    #print contigs according to the input order
    if reference != "":
        print reference.strip()
    for entry in header["CONTIGS"]:
        print(entry.strip())
    del header["CONTIGS"]
    for entry in sorted(header["FILTER"]):
        print(header["FILTER"][entry].strip())
    del header["FILTER"]

    all_formats=[]
    for entry in sorted(header["FORMAT"]):
        print(header["FORMAT"][entry].strip())
        

    del header["FORMAT"]

    #print the other lines in lexiographic order
    for key in sorted(header):
        for entry in sorted(header[key]):
            print(header[key][entry].strip())
    #print subheaders
    for entry in sorted(subheader):
        print(subheader[entry].strip())

    print("##INFO=<ID=SVID,Number=1,Type=String,Description=\"The variant IDs of the replaced SV calls\">")
    print columns
    return ()

parser = argparse.ArgumentParser("""MobileAnn - mobile element classification of SV calls""")
parser.add_argument('--sv', type=str, help="a vcf containing sv", required = True)
parser.add_argument('--me', type=str, help="a vcf containing the mobile elements", required = True)
parser.add_argument('--rm', type=str, help="a repeat masker bed file (format:chr<tab>pos<tab>end<tab>repeat)", required = True)
parser.add_argument('-d', type=int,default=150, help="maximum distance between sv call and mobile element/repeat (default=150)")
args = parser.parse_args()

header=[]
variants=[]

#construct and print the header
construct_header(args)

#load the sv calls
sv_pos=[]
sv_chr=[]
sv_lines=[]
sv_id=[]
contig_order=[]
for line in open(args.sv):
    if line[0] == "#":
        continue
                
    chrA, posA, chrB, posB,event_type,INFO,format = readVCF.readVCFLine(line)

    if not chrA in contig_order:
        contig_order.append(chrA)

    sv_pos.append([posA,posB])
    sv_chr.append([chrA,chrB])
    sv_lines.append(line.strip())
    sv_id.append(line.split()[2])

sv_pos=numpy.array(sv_pos)

#load the repats

repeats={}
first=True
for line in open(args.rm):
    if line[0] == "#":
        continue
    if first:
        first=False
        continue
    content=line.strip().split()
    if not content[0] in repeats:
        repeats[content[0]] = []
    repeats[content[0]].append([int(content[1])-args.d,int(content[2])+args.d])

for chromosome in repeats:
    repeats[chromosome]=numpy.array(repeats[chromosome])

me_lines=[]
conn = sqlite3.connect(":memory:")
c=conn.cursor()
c.execute("CREATE TABLE ME (chr TEXT, pos INT, idx INT)")
lines=[]
i=0

#load the me file
me_lines=[]
for line in open(args.me):
    if line[0] == "#":
        continue
    
    content=line.strip().split()
    chr=content[0]
    pos=int(content[1])
    lines.append([chr,pos,i])
    me_lines.append(line.strip())
    i+=1

c.executemany('INSERT INTO ME VALUES (?,?,?)',lines)
c.execute("CREATE INDEX select_pos on ME (chr,pos,pos)")
conn.commit()
del lines

skip=[]
me_to_print=[]
me_index=[]

#match the repeats MEIs and sv calls
for i in range(0,len(sv_lines)):
    found=False
    repeat=[]
    if sv_chr[i][1] in repeats:
        repeat=repeats[sv_chr[i][1]][numpy.where( (repeats[sv_chr[i][1]][:,0] < sv_pos[i][1]) & (repeats[sv_chr[i][1]][:,1] > sv_pos[i][1]) )]
    if len(repeat):
        A='SELECT idx FROM ME WHERE chr == \'{}\' AND pos < {} AND pos > {} '.format(sv_chr[i][0],sv_pos[i][0]+args.d,sv_pos[i][0]-args.d)
            
        for hit in c.execute(A):
            if not hit[0] in me_index:
                me_index.append(hit[0])
                sv_line=sv_lines[i]
                content=me_lines[int(hit[0])].split()
                content[7]=merge_info(sv_lines[i].split("\t")[7],content[7])
                content[7]+=";SVID={}".format(sv_id[i])
                del content[8:]
                content += sv_line.split("\t")[8:]
                me_to_print.append("\t".join(content))
            found=True

    repeat=[]
    if sv_chr[i][0] in repeats:
        repeat=repeats[sv_chr[i][0]][numpy.where( (repeats[sv_chr[i][0]][:,0] < sv_pos[i][1]) & (repeats[sv_chr[i][0]][:,1] > sv_pos[i][1]) )]
    
    if len(repeat):
        A='SELECT idx FROM ME WHERE chr == \'{}\' AND pos < {} AND pos > {} '.format(sv_chr[i][1],sv_pos[i][1]+args.d,sv_pos[i][1]-args.d)            
        for hit in c.execute(A):
            if not hit[0] in me_index:
                me_index.append(hit[0])
                sv_line=sv_lines[i]
                content=me_lines[int(hit[0])].split()
                content[7]=merge_info(sv_lines[i].split("\t")[7],content[7])
                content[7]+=";SVID={}".format(sv_id[i])
                del content[8:]
                content += sv_line.split("\t")[8:]
                me_to_print.append("\t".join(content))
            
            found=True

    if found:
        skip.append(i)

#order the calls
variants={}
for me in me_to_print:
    content=me.split()
    if not content[0] in variants:
        variants[content[0]]=[]
    content[1]=int(content[1])
    variants[content[0]].append(content)

skip=set(skip)
for i in range(0,len(sv_lines)):
    if not i in skip:
        content=sv_lines[i].split()
        if not content[0] in variants:
            variants[content[0]]=[]
    content[1]=int(content[1])
    variants[content[0]].append(content)

#print sorted calls
for chromosome in contig_order:
    for variant in sorted(variants[chromosome], key=lambda x: x[1]):
        variant[1]=str(variant[1])
        print "\t".join(variant)

