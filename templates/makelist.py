#!/usr/bin/env python

vcfs = "$vcf".split(" ")
with open("vcfs.txt", "w") as fh:
    for vcf in vcfs:
        print(vcf, file=fh)
