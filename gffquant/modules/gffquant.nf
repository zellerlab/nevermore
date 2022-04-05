process run_gffquant {
	publishDir "${params.output_dir}", mode: params.publish_mode

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
	gffquant ${db} ${bam} -o ${sample}/${sample} ${params.gffquant_params} > logs/${sample}.o 2> logs/${sample}.e
	"""
}

process collate_feature_counts {
	publishDir "${params.output_dir}", mode: params.publish_mode

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
