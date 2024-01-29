process fastqc {
    
    input:
    tuple val(sample), path(reads)
    val(stage)

    output:
    tuple val(sample), path("stats/${stage}/fastqc/*fastqc_data.txt"), emit: stats
    tuple val(sample), path("stats/${stage}/read_counts/*.${stage}.txt"), emit: counts

    script:
    def compression = (reads[0].name.endsWith(".gz")) ? "gz" : "bz2"
        
    def fastqc_cmd = "fastqc -t ${task.cpus} --extract --outdir=fastqc"


    r1_files = reads.findAll( { it.name.endsWith("_R1.fastq.${compression}") } )
    r2_files = reads.findAll( { it.name.endsWith("_R2.fastq.${compression}") } )
    
    def process_r1 = ""
    def process_r2 = ""
    def r1_input = ""
    if (r1_files.size() != 0) {
        r1_input += "${r1_files.join(' ')}"
        process_r1 = "${fastqc_cmd} ${r1_input} && mv fastqc/*_R1_fastqc/fastqc_data.txt fastqc/${sample.id}_R1_fastqc_data.txt"
    }
    def r2_input = ""
    if (r2_files.size() != 0) {
        r2_input += "${r2_files.join(' ')}"
        process_r2 = "${fastqc_cmd} ${r2_input} && mv fastqc/*_R2_fastqc/fastqc_data.txt fastqc/${sample.id}_R2_fastqc_data.txt"
    }

    


    // ${fastqc_cmd} ${r1_input} && mv fastqc/*_R1_fastqc/fastqc_data.txt fastqc/${sample.id}_R1_fastqc/${sample.id}_R1_fastqc_data.txt
    """
    set -e -o pipefail
    mkdir -p stats/${stage}/read_counts fastqc/
    ${process_r1}
    ${process_r2}

    grep "Total Sequences" fastqc/*data.txt > seqcount.txt
    echo \$(wc -l seqcount.txt)\$'\t'\$(head -n 1 seqcount.txt | cut -f 2) > stats/${stage}/read_counts/${sample.id}.${stage}.txt
	mv fastqc stats/${stage}/
    """
}
