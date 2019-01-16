#!/bin/bash

# make_mutation_types.sh
# AUTHOR: Yulia Rubanova
# Modified for package TrackSig by Cait Harrigan

if [ ! -f $mutation_types_file ] || [ ! -s  $mutation_types_file ]; then
	echo "Type file..."
	perl $get_mutation_type_script muse $vcf_file $phi_file  $mut_order_file >> $mutation_types_file #2>> $log_dir/log.txt
	rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
fi

if [ ! -f $mutation_types_file ] || [ ! -s  $mutation_types_file ]; then
	echo "ERROR: $mutation_types_file file not created. Aborting..."
	exit -1
fi
