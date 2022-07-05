#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { classify_sample } from "./nevermore/modules/functions"
include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"


workflow {

	fastq_ch = Channel
		.fromPath(params.input_dir + "/" + "**.{fastq.gz,fq.gz}")
		.map { file ->
				def sample = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "")
				sample = sample.replaceAll(/_R?[12]$/, "")
				return tuple(sample, file)
		}
		.groupTuple(sort: true)
        .map { classify_sample(it[0], it[1]) }

	nevermore_main(fastq_ch)

	if (!params.skip_profiling) {

		gffquant_flow(nevermore_main.out.alignments)

	}

}
