#!/bin/bash

# header.sh
# AUTHOR: Yulia Rubanova
# Modified for package TrackSig by Cait Harrigan

packagePath=$1
vcf_file=$2 # Path to directory with vcf files
phi_file=$3 # Path to directory with phi.txt files produces by cals_ssm_phis.py

mutation_counts_path="data/counts/" # Path where to write the mutation counts and sum of quadratic phis
mut_order_path="data/mut_order/" # Path where to write mutation ordering (list of mutations sorted by phi). Needed to run get_clusters_at_timepoints.R. get_clusters_at_timepoints.R   composes list of tree node assignments for chunks of 100 mutations (prevalent tree node assignments at this time point)
mutation_types_path="data/mut_types/" # Path where to write files listing mutation type (out of 96 trinucleotide-based types) for each mutation in vcf file sorted by phi
mutation_bootstrap_path="data/bootstrap/" # Path where to write bootstrapped mutation counts

make_hundreds_script="$packagePath/python/make_hundreds.py"
get_mutation_type_script="$packagePath/perl/getMutationTypes.pl"
bootstrap_mutations_script="$packagePath/python/bootstrap_mutations.py"

do_bootstrap=false

N_BOOTSTRAPS=30

log_dir="."

min_number() {
	printf "%s\n" "$@" | sort -g | head -n1
}


if [[ -z $vcf_file ]]; then
	echo "Please provide VCF path"
	exit
fi

if [[ -z $phi_file ]]; then
	echo "Please provide phi path"
	exit
fi

if [[ -z $mutation_counts_path ]]; then
	mutation_counts_path=.
else
   if [ ! -d "$mutation_counts_path" ]; then
	mkdir -p $mutation_counts_path
   fi
fi

if [[ -z $mut_order_path ]]; then
	mut_order_path=.
else
  if [ ! -d "$mut_order_path" ]; then
	mkdir -p $mut_order_path
  fi
fi

if [[ -z $mutation_types_path ]]; then
	mutation_types_path=.
else
   if [ ! -d "$mutation_types_path" ]; then
	mkdir -p $mutation_types_path
  fi
fi

if [[ -z $mutation_bootstrap_path ]]; then
	mutation_bootstrap_path=$mutation_counts_path
else
   if [ ! -d "$mutation_bootstrap_path" ]; then
	mkdir -p $mutation_bootstrap_path
  fi
fi


if [[ ! -f $vcf_file || $vcf_file != *.vcf ]]; then
	echo "Skipping $vcf_file"
	continue
fi

tumor_id=$(basename ${vcf_file})
tumor_id=${tumor_id%.vcf}

if [[ ! -f "$phi_file" ]]; then
	echo "$phi_file not found. Aborting..."
exit
fi

mutation_counts_file=$mutation_counts_path/$tumor_id.${tumor_part}phi.txt
mutation_types_file=$mutation_types_path/$tumor_id.${tumor_part}mut_types.txt
mutation_quadraticp_file=$mutation_counts_path/$tumor_id.${tumor_part}quadraticp.txt
mut_order_file=$mut_order_path/$tumor_id.${tumor_part}mut_order.txt

