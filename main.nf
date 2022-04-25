#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_simple_preprocessing } from "./nevermore/workflows/nevermore"
include { classify_sample } from "./nevermore/modules/functions"
include { remove_host_kraken2_individual; remove_host_kraken2 } from "./nevermore/modules/decon/kraken2"
include { prepare_fastqs } from "./nevermore/modules/converters/prepare_fastqs"
include { gffquant_flow } from "./gffquant/workflows/gffquant"
include { fastqc } from "./nevermore/modules/qc/fastqc"
include { multiqc } from "./nevermore/modules/qc/multiqc"

def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)

def config_dir = (projectDir.endsWith("nevermore")) ? "${projectDir}/config" : "${projectDir}/nevermore/config"


process collate_stats {
	publishDir params.output_dir, mode: params.publish_mode

	input:
	path(stats_files)

	output:
	path("reports/read_count_table.txt")

	script:
	"""
	mkdir -p reports/
	collate_stats.py . > reports/read_count_table.txt
	"""
}


process bwa_mem_align {
	label 'align'

	input:
	tuple val(sample), path(reads)
	path(reference)

	output:
	tuple val(sample), path("${sample.id}.bam"), emit: bam

	script:
	def align_cpus = 4 // figure out the groovy division garbage later (task.cpus >= 8) ?
	def sort_cpus = 4
	def reads2 = (sample.is_paired) ? "${sample.id}_R2.fastq.gz" : ""
	def blocksize = "-K 100000000"  // shamelessly taken from NGLess

	"""
	bwa mem -a -t ${align_cpus} ${blocksize} \$(readlink ${reference}) ${sample.id}_R1.fastq.gz ${reads2} | samtools view -F 4 -buSh - | samtools sort -@ ${sort_cpus} -o ${sample.id}.bam
	"""

}


process merge_and_sort {
	label 'samtools'
	publishDir params.output_dir, mode: params.publish_mode

	input:
	tuple val(sample), path(bamfiles)

	output:
	tuple val(sample), path("bam/${sample}.bam"), emit: bam
	tuple val(sample), path("stats/bam/${sample}.flagstats.txt"), emit: flagstats

	script:
	// need a better detection for this
	if (bamfiles instanceof Collection && bamfiles.size() >= 2) {
		"""
		mkdir -p bam/ stats/bam
		samtools merge -@ $task.cpus bam/${sample}.bam ${bamfiles}
		samtools flagstats bam/${sample}.bam > stats/bam/${sample}.flagstats.txt
		"""
	} else {
		// i don't like this solution
		"""
		mkdir -p bam/ stats/bam
		ln -s ../${bamfiles[0]} bam/${sample}.bam
		samtools flagstats bam/${sample}.bam > stats/bam/${sample}.flagstats.txt
		"""
	}
}


process merge_single_fastqs {
	input:
	tuple val(sample), path(fastqs)

	output:
	tuple val(sample), path("merged/${sample.id}_R1.fastq.gz"), emit: fastq

	script:
	"""
	mkdir -p merged/

	cat *.fastq.gz > merged/${sample.id}_R1.fastq.gz
	"""

}


workflow nevermore_align {

	take:

		fastq_ch

	main:

		/*	route all single-read files into a common channel */

		single_ch = fastq_ch
			.filter { it[0].is_paired == false }
			.map { sample, fastq ->
				def meta = [:]
				meta.id = fastq.name.replaceAll(/_R1.fastq.gz$/, "")
				meta.is_paired = false
				meta.merged = false
				return tuple(meta, fastq)
			}

		/*	route all paired-end read files into a common channel */

		paired_ch = fastq_ch
			.filter { it[0].is_paired == true }
			.map { sample, fastq ->
				def meta = [:]
				meta.id = sample.id
				meta.is_paired = true
				meta.merged = true
				return tuple(meta, fastq)
			}

		/*	group all single-read files by sample and route into merge-channel */

		merged_single_ch = single_ch
			.map { sample, fastq ->
				return tuple(
					sample.id.replaceAll(/.(orphans|singles|chimeras)$/, ".singles"),
					fastq
				)
			}
			.groupTuple(sort: true)
			.map { sample_id, files ->
				def meta = [:]
				meta.id = sample_id
				meta.is_paired = false
				meta.merged = true
				return tuple(meta, files)
			}

		/*	then merge single-read file groups into single files */

		merge_single_fastqs(merged_single_ch)

		/* 	take all single-read files except for the qc-survivors,
			concat with merged single-read files (takes care of single-end qc-survivors),
			concat with paired-end files,
			and route them into a channel for post-qc fastqc analysis
		*/

		fastqc_in_ch = single_ch
			.filter { ! it[0].id.endsWith(".singles") }
			.map { sample, fastq ->
				def meta = [:]
				meta.id = fastq.name.replaceAll(/_R1.fastq.gz$/, "")
				meta.is_paired = false
				meta.merged = false
				return tuple(meta, fastq)
			}
			.concat(merge_single_fastqs.out.fastq)
			.concat(paired_ch)

		/*	perform post-qc fastqc analysis and generate multiqc report on merged single-read and paired-end sets */

		fastqc(fastqc_in_ch, "qc")

		multiqc(
			fastqc.out.stats
				.filter { it[0].merged == true || it[0].is_paired == true }
				.map { sample, report -> return report }
				.collect(),
			"${config_dir}/multiqc.config",
			"qc"
		)

		/*	align merged single-read and paired-end sets against reference */

		bwa_mem_align(
			paired_ch.concat(merge_single_fastqs.out.fastq),
			params.reference
		)

		/*	merge paired-end and single-read alignments into single per-sample bamfiles */

		aligned_ch = bwa_mem_align.out.bam
			.map { sample, bam ->
				sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_id, bam)
			}
			.groupTuple(sort: true)

		merge_and_sort(aligned_ch)


	emit:

		alignments = merge_and_sort.out.bam
		read_counts = fastqc.out.counts
		aln_counts = merge_and_sort.out.flagstats

}


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


	if (do_preprocessing) {

		prepare_fastqs(fastq_ch)

		raw_fastq_ch = prepare_fastqs.out.reads

		nevermore_simple_preprocessing(raw_fastq_ch)

		preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out
		if (!params.drop_orphans) {
			preprocessed_ch = preprocessed_ch.concat(nevermore_simple_preprocessing.out.orphan_reads_out)
		}

		if (params.remove_host) {

			remove_host_kraken2_individual(preprocessed_ch, params.remove_host_kraken2_db)

			preprocessed_ch = remove_host_kraken2_individual.out.reads
			if (!params.drop_chimeras) {
				chimera_ch = remove_host_kraken2_individual.out.chimera_orphans
					.map { sample, file ->
						def meta = [:]
						meta.is_paired = false
						meta.id = sample.id + ".chimeras"
						return tuple(meta, file)
					}
				preprocessed_ch = preprocessed_ch.concat(chimera_ch)
			}

		}

	} else {

		preprocessed_ch = fastq_ch

	}


	nevermore_align(preprocessed_ch)

	gffquant_flow(nevermore_align.out.alignments)

	if (do_preprocessing) {

		collate_ch = nevermore_simple_preprocessing.out.read_counts
			.map { sample, file -> return file }
			.collect()
			.concat(
				nevermore_align.out.read_counts
					.map { sample, file -> return file }
					.collect()
			)
			.concat(
				nevermore_align.out.aln_counts
					.map { sample, file -> return file }
					.collect()
			)
			.collect()
		// collate_ch.view()
		collate_stats(
			collate_ch
		)

	}
}
