## Order of RRBS analysis pipeline<br>
#Order	scripts<br>
1	<td>run_fastqc.sh (batch)<br>
2	generate_scripts.sh & [sample]_trim.sh (batch)<br>
3	run_bwa_meth_indexing.sh<br>
4	generate_scripts.sh & [sample]_align.sh (batch)<br>
5	bismark_genome_index.sh<br>
6	generate_scripts.sh & [sample]_bismark.sh (batch)<br>
7	make_methylKit_objects.R & make_methylKit_objects.sh (batch)<br>
8	methyKit_DMRs_new_final_metadata.R & run_methyKit_DMRs.sh (batch)<br>
