#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bam_input } from "./nevermore/workflows/input"
include { run_gffquant: collate_feature_counts } from "./nevermore/workflows/gffquant"

if (params.input_dir && params.remote_input_dir) {
	log.info """
		Cannot process both --input_dir and --remote_input_dir. Please check input parameters.
	""".stripIndent()
	exit 1
} else if (!params.input_dir && !params.remote_input_dir) {
	log.info """
		Neither --input_dir nor --remote_input_dir set.
	""".stripIndent()
	exit 1
}

def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir

workflow {

	if (!params.bam_input_pattern) {
		bam_input(
			Channel.fromPath(input_dir + "/**.bam", type: "file")
		)
	} else if (bam_input_pattern) {
		bam_input_pattern = input_dir + "/" + params.bam_input_pattern
		bam_input(
			Channel.fromPath(bam_input_pattern)
		)
	}

	alignment_ch = bam_input.out.bamfiles
		.map { sample, bam ->
			return tuple(sample.id, bam[0])
		}
		.groupTuple(sort: true)

	run_gffquant(
		alignment_ch,
		params.gffquant_db
	)

	feature_count_ch = run_gffquant.out.results
		.map { sample, files -> return files }
		.flatten()
		.filter { !it.name.endsWith("Counter.txt.gz") }
		.filter { params.collate_gene_counts || !it.name.endsWith("gene_counts.txt.gz") }
		.map { file ->
			def category = file.name
				.replaceAll(/\.txt\.gz$/, "")
				.replaceAll(/.+\./, "")
			return tuple(category, file)
		}
		.groupTuple(sort: true)
		.combine(
			Channel.from(params.gq_collate_columns.split(","))
		)

	feature_count_ch.view()

	if (!params.no_collate) {
		collate_feature_counts(feature_count_ch)
	}

}
