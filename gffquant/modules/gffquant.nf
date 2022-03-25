process run_gffquant {
	publishDir "${output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(bam)
	path(db)

	output:
	tuple val(sample), path("${sample}/*.txt.gz"), emit: results

	script:
	def emapper_version = (params.emapper_version) ? "--emapper_version ${params.emapper_version}" : ""
	"""
	echo $sample $bam
	mkdir -p logs/
	gffquant ${db} ${bam} -o ${sample}/${sample} -m ${params.mode} --ambig_mode ${params.ambig_mode} ${emapper_version} ${params.strand_specific} > logs/${sample}.o 2> logs/${sample}.e
	"""
}

process collate_feature_counts {
	publishDir "${output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(count_tables)

	output:
	path("collated/*.txt.gz"), emit: collated, optional: true

	script:
	"""
	mkdir -p collated/
	collate_counts . -o collated/collated -c uniq_scaled
	collate_counts . -o collated/collated -c combined_scaled
	"""
}


workflow gffquant_flow {

	take:

		bam_ch

	main:

		run_gffquant(bam_ch, params.gffquant_db)

		feature_count_ch = run_gffquant.out.results //.collect()
			.map { sample, files -> return files }
			.flatten()
			.filter { !it.name.endsWith("gene_counts.txt") }
			.filter { !it.name.endsWith("seqname.uniq.txt") }
			.filter { !it.name.endsWith("seqname.dist1.txt") }
			.map { file -> 
				def category = file.name.replaceAll(/\.txt$/, "")
					.replaceAll(/.+\./, "")
				return tuple(category, file)
			}
			.groupTuple(sort:true)

		//feature_count_ch.view()
	
		collate_feature_counts(feature_count_ch)

	emit:

		counts = run_gffquant.out.results
		collated = collate_feature_counts.out.collated

}