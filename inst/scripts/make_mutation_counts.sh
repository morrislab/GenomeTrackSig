#!/bin/bash

# make_mutation_counts.sh
# AUTHOR: Yulia Rubanova
# Modified for package TrackSig by Cait Harrigan

num_mutations=$(min_number `cat $phi_file | wc -l` `cat $vcf_file | wc -l`)
num_hundreds=$(($num_mutations/100 + ($num_mutations % 100 > 0)))


if [ ! -f $mutation_counts_file ]; then
	if [ $num_mutations -lt 100 ]; then
		echo "Less than 100 mutaions in a file $vcf_file or $phi_file"
		touch $mutation_counts_file
	else
		echo "Count file..."
			for i in `seq 1 $num_hundreds`; do
			if [ $num_mutations -ge $((i*100-1)) ]; then
					###echo $i $((i*100-1))
				python $make_hundreds_script $mutation_types_file  $((i*100-100)) $((i*100-1)) |\
				awk '{split($0, line, ";"); print line[1] >> out1; print line[2] >> out2}'  out1="$mutation_counts_file" out2="$mutation_quadraticp_file"
				rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
				fi
		done
	fi
fi


if [[ ! -f "$mutation_counts_file" ]]; then
	echo "ERROR: $mutation_counts_file file not created. Aborting..."
	exit -1
fi

if [ "$do_bootstrap" = true ] ; then
   bootstrap_dir=$mutation_bootstrap_path/$tumor_id/
   if [ !  -z  $tumor_part  ]; then
	bootstrap_dir=${bootstrap_dir:0:${#bootstrap_dir}-1}
	bootstrap_dir=$bootstrap_dir.${tumor_part:0:${#tumor_part}-1}/
   fi

   if [ ! -d $bootstrap_dir ]; then
	mkdir $bootstrap_dir
   fi

   for i in `seq 1 $N_BOOTSTRAPS`; do
	mutation_bootstrap_file=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut.txt
	mutation_bootstrap_counts_file=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut_counts.txt
	mutation_bootstrap_quadraticp_file=$bootstrap_dir/$tumor_id.${tumor_part}$i.quadraticp.txt

	mutation_bootstrap_file_unsorted=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut.unsorted.txt
	mutation_bootstrap_counts_file_unsorted=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut_counts.unsorted.txt

	if [ ! -f $mutation_bootstrap_file_unsorted ] || [ ! -s  $mutation_bootstrap_file_unsorted ] || [ ! -f $mutation_bootstrap_file ] || [ ! -s  $mutation_bootstrap_file ]; then

		if [ $num_mutations -lt 100 ]; then
					touch $mutation_bootstrap_file
			else
			    # echo "Bootstrap mutations..."
				python $bootstrap_mutations_script $mutation_types_file  $num_mutations  > $mutation_bootstrap_file_unsorted #2>>$log_dir/log.txt
				rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

				sort -k3 -nr $mutation_bootstrap_file_unsorted > $mutation_bootstrap_file
		fi
	fi

	 if [ ! -f $mutation_bootstrap_counts_file ] || [ ! -s  $mutation_bootstrap_counts_file ]; then
			# echo "Bootstrap counts..."
			for t in `seq 1 $num_hundreds`; do
				if [ $num_mutations -ge $((t*100-1)) ]; then
					python $make_hundreds_script $mutation_types_file  $((t*100-100)) $((t*100-1)) |\
					awk '{split($0, line, ";"); print line[1] >> out1; print line[2] >> out2}'  out1="$mutation_bootstrap_counts_file" out2="$mutation_bootstrap_quadraticp_file"
					rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
				fi
			done
		fi

	# if [ ! -f $mutation_bootstrap_counts_file_unsorted ] || [ ! -s  $mutation_bootstrap_counts_file_unsorted ]; then
	# 		#echo "Bootstrap counts on unsorted mutations..."
	# 		for t in `seq 1 $num_hundreds`; do
	# 			if [ $num_mutations -ge $((t*100-1)) ]; then
	# 				python $make_hundreds_script $mutation_bootstrap_file_unsorted  $((t*100-100)) $((t*100-1)) >> $mutation_bootstrap_counts_file_unsorted # 2>>$log_dir/log.txt
	# 			fi
	# 		done
	# fi

  done
fi

