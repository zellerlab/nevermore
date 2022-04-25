import argparse
import os
import re
import sys


def collect_readcounts(d):
	readcounts = {}
	if os.path.exists(os.path.join(d, "raw_counts")):
		walk = os.walk(os.path.join(d, "raw_counts"))
		pwd, dirs, files = next(walk)

		for f in files:
			prefix, suffix = os.path.splitext(f)
			if suffix == ".txt":
				prefix, suffix = os.path.splitext(prefix)
				if suffix not in (".chimeras", ".singles", ".orphans"):
					prefix, suffix = prefix + suffix, ".main"

				count_data = open(os.path.join(pwd, f)).read().strip()
				n_mates = int(count_data[0])
				n_reads = int(count_data.split("\t")[1])

				readcounts.setdefault(prefix, {})[suffix] = n_mates * n_reads
				if suffix == ".main":
					readcounts[prefix]["paired_end"] = n_mates == 2
					readcounts[prefix]["single_end"] = readcounts[prefix].get("single_end", False) or n_mates == 1
				if suffix == ".singles":
					readcounts[prefix]["single_end"] = True
				

	return readcounts


def collect_flagstats(d):

	def extract(line):
		line = line.strip().split(" ")
		return int(line[0]) + int(line[2])

	flagstats = {}
	if os.path.exists(os.path.join(d, "stats")):
		walk = os.walk(os.path.join(d, "stats"))
		pwd, dirs, files = next(walk)

		for f in files:
			lines = open(os.path.join(pwd, f)).readlines()
			sample = f.replace(".flagstats.txt", "")
			flagstats[sample] = {
				"primary_alignments": extract(lines[0x1]),
				"secondary_alignments": extract(lines[0x2]),
				"npaired_mapped": extract(lines[0xC]),
				"nsingles_mapped": extract(lines[0xD]) + (extract(lines[0x7]) - extract(lines[0x8])),
			}
			flagstats[sample]["n_alignments"] = flagstats[sample]["primary_alignments"] + flagstats[sample]["secondary_alignments"]

	return flagstats


def write_output(readcounts, flagstats, outfile):
	samples = sorted(set(readcounts).union(flagstats))

	header = ["sample", "pe", "se", "n_reads", "main_lib", "add_lib", "qc_orphans", "dc_orphans"]
	header += ["n_aln", "pri_aln", "sec_aln", "pe_aln", "se_aln"]
	print(*header, sep="\t", file=outfile)

	for sample in samples:

		if readcounts[sample].get(".main") is None:
			readcounts[sample][".main"] = readcounts[sample][".singles"]
			del readcounts[sample][".singles"]


		line = [sample, readcounts[sample].get("paired_end", False), readcounts[sample].get("single_end", False), sum(v for k, v in readcounts[sample].items() if k[0] == ".")]
		line += (readcounts[sample].get(key, 0) for key in (".main", ".singles", ".orphans", ".chimeras"))
		line += (flagstats[sample].get(key, 0) for key in ("n_alignments", "primary_alignments", "secondary_alignments", "npaired_mapped", "nsingles_mapped"))

		print(*line, sep="\t", file=outfile)



def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("input_dir", type=str, default=".")
	args = ap.parse_args()

	readcounts = collect_readcounts(args.input_dir)
	flagstats = collect_flagstats(args.input_dir)

	with open(os.path.join(args.input_dir, "reports", "library_stats.tsv"), "wt") as _out:
		write_output(readcounts, flagstats, _out)
	

if __name__ == "__main__":
	main()

"""
0	299821609 + 0 in total (QC-passed reads + QC-failed reads)
1	51499902 + 0 primary
2	244967945 + 0 secondary
3	3353762 + 0 supplementary
4	0 + 0 duplicates
5	0 + 0 primary duplicates
6	299821609 + 0 mapped (100.00% : N/A)
7	51499902 + 0 primary mapped (100.00% : N/A)
8	47968663 + 0 paired in sequencing
9	23987824 + 0 read1
A	23980839 + 0 read2
B	23298984 + 0 properly paired (48.57% : N/A)
C	46949898 + 0 with itself and mate mapped
D	1018765 + 0 singletons (2.12% : N/A)
E	23643000 + 0 with mate mapped to a different chr
F	9277385 + 0 with mate mapped to a different chr (mapQ>=5)
"""
