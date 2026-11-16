#!/usr/bin/env python
#Python script to filter vcf file based on allelic ratio of coverage
#Fitlering is done on an individual bases
#Genotype not passing the filter are assigned missing value
#usage: filtering_vcf.py --vcf <file> --lratio <float> -- hratio <float>

import gzip as gz
import argparse as arg
from collections import defaultdict

## Define class
class vcf_parser_line():
    sample_list = []
    #filtering_status = defaultdict(lambda: defaultdict(int))

    def __init__(self, vcf_line: str):
        self.line = vcf_line.strip().split("\t")

        if vcf_line.startswith("#CHROM"):
            vcf_parser_line.sample_list = self.line[9:]

        if vcf_line.startswith("#"):
            self.header = True
            return

        self.header = False
        self.linestart = self.line[0:8]
        self.info = self.line[8].split(":")
        self.info_dict = defaultdict(list)

        for ind in self.line[9:]:
             for pos, field in enumerate(ind.split(":")):
                 self.info_dict[self.info[pos]].append(field)


    def rebuild_line(self):
        geno = "\t".join([":".join(x) for x in zip(*self.info_dict.values())])
        return "\t".join(self.linestart) + "\t" + ":".join(self.info) + "\t" + geno

    def filter_ratio(self, thresh):
        if "AD" not in self.info:
            print("cannot filter on ratio")
            return

        for ind in range(len(vcf_parser_line.sample_list)):
            AD = [int(x) for x in self.info_dict["AD"][ind].split(",")]
            AD = [x / sum(AD) if sum(AD) >0 else 0 for x in AD]

            if (max(AD) <= thresh[1] and min(AD) >= thresh[0]):
                self.info_dict["GT"][ind] = "."


def myopen(file, mode = "rt"):
    if file.endswith(".gz"):
        return gz.open(file, mode)
    return open(file, mode)


## Define argument parser
parser = arg.ArgumentParser()
parser.add_argument("--vcf", type = str, required = True)
parser.add_argument("--lratio", default = 0, type = float)
parser.add_argument("--hratio", default = 1, type = float)
args = parser.parse_args()

out = "06_vcf/ratio_qual_filtered" + "_" + str(args.lratio)+ "_" + str(args.hratio) + ".vcf"
indat = myopen(args.vcf)
output = myopen(out, "wt")

for line in indat:
    l = vcf_parser_line(line)
    if l.header:
        output.write(line)
        continue

    if l.linestart[4] != "." : #check if polymorpic, then filter genotype
        l.filter_ratio([args.lratio, args.hratio])
    output.write(l.rebuild_line() + "\n")

