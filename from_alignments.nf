#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bam_input } from "./nevermore/workflows/input"
include { gffquant_flow } from "./nevermore/workflows/gffquant"

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

	gffquant_flow(
		alignment_ch,
	)

}
