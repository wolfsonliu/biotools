#! /bin/env python3


name_index = {}

with open("knowgene", "r") as index:
    for line in index:
        if not line.startswith("#"):
            name_index[line.split("\t")[0]] = line.split("\t")[-1].rstrip()

name_change = []

with open("result.tsv", "r") as destfile:
    name_change = [line.rstrip() for line in destfile]

output = open("result_gs.tsv", "w")

for record in name_change:
    if record.startswith("id"):
        output.write(record + "\tgenesymbol\n")
    else:
        output.write(record + "\t" + name_index[record.split("\t")[0]] + "\n")

output.close()
