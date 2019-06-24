# AUTHOR: Yulia Rubanova

#!/usr/bin/python

from __future__ import print_function
import sys
import pandas as pd

NUCLEOTIDES = ['A', 'C', 'T', 'G']

def make_bins(vcaf, start, end):

  print(vcaf.type())
	file = sys.argv[1]
	start = sys.argv[2]
	if start.isdigit():
		start = int(start)

	end = sys.argv[3]
	if end.isdigit():
		end = int(end)

	assert (end >= start)

	counts = {}
	for first in NUCLEOTIDES:
		for third in NUCLEOTIDES:
			for ref in ['C', 'T']:
				for alt in NUCLEOTIDES:
					if ref != alt:
						key = ref + "_" + alt + "_" + first + ref + third
						counts[key] = 0

	sum_phis = 0
	sum_squ_phis = 0

	with open(file, 'r') as mut_types:
		for i, line in enumerate(mut_types):
			if i < start or i > end:
				continue

			try:
				chr, pos, phi, ref, alt, context = line.split()
			except:
				print("Couldn't parse the line: " + line)
				continue
			if  phi == "NA":
				continue

			counts[ref + "_" + alt + "_" + context] += 1
			sum_phis += float(phi)
			sum_squ_phis += float(phi) ** 2

	assert(sum(counts.values()) == max(0,end - start + 1))

	mean_phi = sum_phis / (end - start + 1)
	mean_squ_phi = sum_squ_phis / (end - start + 1)

	print(file + "\t", end="")
	print(str(mean_phi) + "\t", end="")
	for mutation in sorted(counts.keys()):
		print(str(counts[mutation]) + "\t", end="")

	#append sum square phis to output - use delim ";"
	print(";", end="")
	print(mean_squ_phi, end="")

	print("")





