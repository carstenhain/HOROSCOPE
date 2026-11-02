#!/usr/bin/env nextflow

params.samplesheet = "samplesheet.csv"
params.outdir = null
params.kmer_fasta = "${projectDir}/data/all_tagging_kmers.fasta.gz"
params.kmer_info_file = "${projectDir}/data/all_tagging_kmers.tsv.gz"
params.model_directory = "${projectDir}/models/"
params.help = false

// Help message
def helpMessage() {
    log.info"""
    Usage:
        nextflow run horoscope.nf --samplesheet <path> --outdir <path>

    Required arguments:
        --samplesheet PATH    Path to the semicolon-separated samplesheet file
        --outdir PATH         Directory where result files are stored

    Optional arguments:
        --help                Show this help message and exit

    The samplesheet should be a semicolon-separated file with columns:
      - NAME:                 Sample identifier (string)
      - FILE_PATH:            Path to the input file (string)
      - FILE_TYPE:            Type of the input file (string, e.g. 'bam', 'cram', 'fastq', 'dbg')
      - CRAM_REFERENCE_PATH:  Path to reference FASTA for CRAM decoding (or NA)
      - MAKE_DBG:             Whether to build a de Bruijn graph for this sample (true/false)
      - NORMALIZATION:        Normalization mode ('global', 'chromosome', 'p_arm', 'q_arm', or a float value)
    """.stripIndent()
}

// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

process DBG_LIGHTER {
    tag "LIGHTER_${meta.id}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.id}.cor.fq.gz")

    script:
    """
    lighter_input="${input_file}"

    if [[ "${meta.file_type}" == "bam" ]]; then    
        samtools fastq -@ ${task.cpus} "${input_file}" > "${meta.id}.for_lighter.fastq"
        lighter_input="${meta.id}.for_lighter.fastq"
    elif [[ "${meta.file_type}" == "cram" ]]; then
        if [[ "${meta.cram_ref}" == "NA" || -z "${meta.cram_ref}" ]]; then
            echo "CRAM input requires CRAM_REFERENCE_PATH for sample ${meta.id}" >&2
            exit 1
        fi
        samtools fastq -@ ${task.cpus} --reference "${meta.cram_ref}" "${input_file}" > "${meta.id}.for_lighter.fastq"
        lighter_input="${meta.id}.for_lighter.fastq"
    elif [[ "${meta.file_type}" != "fastq.gz" && "${meta.file_type}" != "fastq" && "${meta.file_type}" != "fasta" ]]; then
        echo "Unsupported FILE_TYPE for DBG_LIGHTER: ${meta.file_type}. Allowed: fastq.gz, fastq, fasta, bam, cram" >&2
        exit 1
    fi

    # store lighter output in this directory
    mkdir -p lighter_output

    # run lighter
    lighter -t ${task.cpus} -r \${lighter_input} -od lighter_output -trim -discard -k 23 3100000000 0.188

    # Capture the expected corrected filename from lighter input naming convention.
    source_name=\$(basename "\${lighter_input}")
    corrected_file=""
    case "\${source_name}" in
        *.fastq.gz)
            corrected_file="lighter_output/\${source_name%.fastq.gz}.cor.fq.gz"
            ;;
        *.fq.gz)
            corrected_file="lighter_output/\${source_name%.fq.gz}.cor.fq.gz"
            ;;
        *.fastq)
            corrected_file="lighter_output/\${source_name%.fastq}.cor.fq"
            ;;
        *.fq)
            corrected_file="lighter_output/\${source_name%.fq}.cor.fq"
            ;;
    esac

    if [[ -z "\${corrected_file}" ]]; then
        echo "Lighter did not produce a corrected read file for sample ${meta.id}" >&2
        exit 1
    fi
    if [[ ! -f "\${corrected_file}" ]]; then
        echo "Corrected read file path was set but file does not exist: \${corrected_file}" >&2
        exit 1
    fi
    if [[ "\${corrected_file}" == *.gz ]]; then
        cp "\${corrected_file}" "${meta.id}.cor.fq.gz"
    else
        gzip -c "\${corrected_file}" > "${meta.id}.cor.fq.gz"
        rm -f "\${corrected_file}"
    fi

    if [[ -f "${meta.id}.for_lighter.fastq" ]]; then
        rm -f "${meta.id}.for_lighter.fastq"
    fi

    """
}

process DBG_BCALM {
    tag "BCALM_${meta.id}"

    input:
    tuple val(meta), path(corrected_fastq)

    output:
    tuple val(meta), path("${meta.id}.unitigs.fa.gz")

    script:
    """
    bcalm \
        -in "${corrected_fastq}" \
        -kmer-size 61 \
        -abundance-min 3 \
        -max-memory ${task.memory.toMega()} \
        -out "${meta.id}"

    gzip -c "${meta.id}.unitigs.fa" > "${meta.id}.unitigs.fa.gz"

    rm -f "${meta.id}.list_reads"
    """
}

process UNPACK_KMER_FASTA {
    input:
    path(kmer_fasta_gz)

    output:
    path("all_tagging_kmers.fasta")

    script:
    """
    gzip -dc "${kmer_fasta_gz}" > all_tagging_kmers.fasta
    """
}

process GATHER_KMERS_DBG {
    tag "${meta.id}"

    input:
    tuple val(meta), path(input_file)
    path(kmer_fasta)
    path(dbg_script)

    output:
    tuple val(meta), path("${meta.id}.kmer.tsv")

    script:
    """
    ### rewrite fasta to a list of kmers
    cat "${kmer_fasta}" | awk 'NR % 2 == 0' > kmer_list.txt

    ### create local uncompressed copy of the dbg file in task work directory
    gzip -dc "${input_file}" > dbg.uncompressed.fa

    ### subset dbg to those contigs with one matching kmers and write to new file
    pv dbg.uncompressed.fa | grep -B 1 -f kmer_list.txt > "${meta.id}".kmers.zgrep.txt
    

    ### call python file to convert the zgrep output into a kmer table
    python3 ${dbg_script} \\
        --name "${meta.id}" \\
        --zgrep_file "${meta.id}.kmers.zgrep.txt"  \\
        --kmer_list kmer_list.txt \\
        --output "${meta.id}.kmer.tsv"

    rm kmer_list.txt
    rm dbg.uncompressed.fa
    """
}

process GATHER_KMERS_JF {
    tag "${meta.id}"

    input:
    tuple val(meta), path(input_file)
    path(kmer_fasta)

    output:
    tuple val(meta), path("${meta.id}.kmer.tsv")

    script:
    """
    if [[ "${meta.file_type}" == "fastq" || "${meta.file_type}" == "fasta" ]]; then
        if [[ "${input_file}" == *.gz || "${input_file}" == *.bgz ]]; then
            jellyfish count -m 61 -s 2G -t ${task.cpus} --if "${kmer_fasta}" -C <(zcat "${input_file}") -o "${meta.id}.jf"
        else
            jellyfish count -m 61 -s 2G -t ${task.cpus} --if "${kmer_fasta}" -C <(cat "${input_file}") -o "${meta.id}.jf"
        fi
    elif [[ "${meta.file_type}" == "bam" ]]; then
        jellyfish count -m 61 -s 2G -t ${task.cpus} --if "${kmer_fasta}" -C <(samtools fasta "${input_file}") -o "${meta.id}.jf"
    elif [[ "${meta.file_type}" == "cram" ]]; then
        if [[ "${meta.cram_ref}" == "NA" || -z "${meta.cram_ref}" ]]; then
            echo "CRAM input requires CRAM_REFERENCE_PATH for sample ${meta.id}" >&2
            exit 1
        fi
        jellyfish count -m 61 -s 2G -t ${task.cpus} --if "${kmer_fasta}" -C <(samtools fasta --reference "${meta.cram_ref}" "${input_file}") -o "${meta.id}.jf"
    else
        echo "Unsupported FILE_TYPE for GATHER_KMERS_JF: ${meta.file_type}" >&2
        exit 1
    fi

    jellyfish query \\
        "${meta.id}.jf" \\
        -s "${kmer_fasta}" \\
        -o "${meta.id}.kmer.tsv"

    rm "${meta.id}.jf"
    """
}

process FORMAT_AND_MERGE_KMERS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(kmer_files)
    path(reformat_script)
    path(kmer_info_file)

    output:
    path("final_kmer_merged.tsv")

    script:
    """
    ls -1 ${kmer_files} > all_kmer_files.txt
    python3 ${reformat_script} \\
        --kmer_file_list all_kmer_files.txt \\
        --kmer_info_file ${kmer_info_file}
    """
}

process GENOTYPE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(merged_kmers)
    path(samplesheet)
    path(genotype_script)

    output:
    path("normalization_metrics.tsv")
    path("centromere_genotyping.tsv")
    path("centromere_genotyping.full.tsv")
    path("centromere_genotyping_results.tsv")

    script:
    """
    python ${genotype_script} \\
        --kmer_table ${merged_kmers} \\
        --samplesheet ${samplesheet} \\
        --model_directory ${params.model_directory}
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    if (params.help) {
        helpMessage()
        exit(0)
    }

    // Parse the semicolon-delimited samplesheet into [meta, input_file] tuples
    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ';')
        .map { row ->
            def normalizationRaw = (row.NORMALIZATION ?: '').toString().trim()
            def normalizationLower = normalizationRaw.toLowerCase()
            def isFloatNormalization = normalizationRaw.isNumber()
            if (!(normalizationLower in ['global', 'chromosome', 'p_arm', 'q_arm']) && !isFloatNormalization) {
                error "Invalid NORMALIZATION value '${normalizationRaw}' for sample ${row.NAME}. Allowed values are 'global', 'chromosome', 'p_arm', 'q_arm', or a float."
            }

            // NORMALIZATION stays out of meta: GENOTYPE reads it from the samplesheet,
            // and keeping it out avoids invalidating cached k-mer counting when it changes
            def meta = [
                id:        row.NAME,
                file_type: row.FILE_TYPE,
                cram_ref:  row.CRAM_REFERENCE_PATH,
                make_dbg:  row.MAKE_DBG.toLowerCase() == 'true'
            ]
            tuple(meta, file(row.FILE_PATH))
        }
        .branch { meta, _input_file ->
            to_build_dbg:  meta.make_dbg
            existing_dbg:  !meta.make_dbg && meta.file_type == 'dbg'
            other:         !meta.make_dbg && meta.file_type != 'dbg'
        }
        .set { branched }

    // Unpack tagging kmers once and reuse in all k-mer gathering tasks
    UNPACK_KMER_FASTA(channel.value(file(params.kmer_fasta, checkIfExists: true)))

    // Build de Bruijn graphs where requested
    DBG_LIGHTER(
        branched.to_build_dbg
    )

    DBG_BCALM(
        DBG_LIGHTER.out
    )

    // Combine freshly built DBGs and pre-existing DBG files, then gather k-mers
    GATHER_KMERS_DBG(
        DBG_BCALM.out.mix(branched.existing_dbg),
        UNPACK_KMER_FASTA.out,
        file("${projectDir}/scripts/dbg_to_kmer_table.py")
    )

    // All remaining file types go to the Jellyfish-based k-mer counter
    GATHER_KMERS_JF(
        branched.other,
        UNPACK_KMER_FASTA.out
    )

    // Collect names of all kmer output files into a single manifest
    FORMAT_AND_MERGE_KMERS(
        GATHER_KMERS_DBG.out.mix(GATHER_KMERS_JF.out)
            .map { _meta, kmer_file -> kmer_file }
            .collect(),
        file("${projectDir}/scripts/reformat_merge_kmers.py"),
        file(params.kmer_info_file, checkIfExists: true)
    )

    GENOTYPE(
        FORMAT_AND_MERGE_KMERS.out,
        channel.value(file(params.samplesheet, checkIfExists: true)),
        file("${projectDir}/scripts/genotype.py")
    )
}
