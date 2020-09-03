params.help = false
if (params.help) {
    log.info """
    -----------------------------------------------------------------------

    freebayes-nf: a better freebayes-parallel
    =========================================

    Output is written to <outdir>/freebayes/<project>.vcf.gz and represents
    a decomposed and normalized VCF and its index.

    required
    --------
    --alignments   Aligned sequences in .bam and/or .cram format. Indexes
                   (.bai/.crai) must be present.
    --fasta        Reference FASTA. Index (.fai) must exist in same
                   directory.

    options
    -------
    --outdir       Base results directory for output. Default: '/.results'
    --project      File prefix for merged and annotated VCF files.
                   Default: 'variants'
    --width        The genomic window size per variant calling job.
                   Default: 1000000
    --options      Arguments to be passed to freebayes command in addition
                   to those already supplied like `--bam`, `--region`, and
                   `--fasta-reference`. Single quote these when specifying
                   on the command line, e.g. --options '--pooled-discrete'.
    --intervals    Picard-style intervals file to use rather than intervals
                   defined in .fai. Something like Broad's interval lists
                   work here if you want to omit masked regions.
    --exclude      Chromosome patterns to omit from variant calling.
                   Default: 'decoy,random,Un,alt,EBV,M,HLA,phi'

    -----------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

// required arguments
params.alignments = false
if( !params.alignments ) { exit 1, "--alignments is not defined" }
params.fasta = false
if( !params.fasta ) { exit 1, "--fasta is not defined" }

// options
params.intervals = false
if( params.intervals ){
    intervals = file(params.intervals)
}
exclude = params.exclude.tokenize(',')
exclude = exclude
    .collect {"$it"}
    .join("|")

// files
fasta = file(params.fasta)
faidx = file("${params.fasta}.fai")

intervals_ch = Channel
    .from(params.intervals ? intervals : faidx)
    .splitCsv(sep: '\t')
    .filter { row ->
        if (row[0] =~ /(${exclude})/) {
            // println("Excluding ${row[0]}")
        } else {
            row
        }
    }
    .map { row ->
        // rows of interval lists
        if (row[0][0] != "@") {
            def interval_start = row[1].toLong()
            def interval_length = row[2].toLong()
            long start
            long end
            int width = params.width

            if (!params.intervals) {
                // update interval start and length for .fai
                interval_start = 0
                interval_length = row[1].toLong()
            }

            while(interval_start < interval_length) {
                start = interval_start
                // add a slight overlap
                end = interval_start + width + 1000
                interval_start = end - 1000
                if (end > interval_length) {
                    end = interval_length
                    interval_start = end
                }
                // add the interval to the channel
                intervals_ch.bind( "${row[0]}:${start}-${end}" )
            }
        }
    }

Channel
    .fromPath(params.alignments, checkIfExists: true)
    .set { alignments_ch }

Channel
    .fromPath(params.alignments, checkIfExists: true)
    .map { file -> file + ("${file}".endsWith('.cram') ? '.crai' : '.bai') }
    .set { indexes_ch }


process run_freebayes {
    tag "$interval"

    input:
    file(aln) from alignments_ch.collect()
    file(idx) from indexes_ch.collect()
    each interval from intervals_ch
    file(fasta)
    file(faidx)

    output:
    // need to rename due to colon in intervals
    file("${params.project}_${interval.replaceAll(~/\:|\-|\*/, "_")}.vcf.gz") into (vcf_ch, makelist_ch)
    file("${params.project}_${interval.replaceAll(~/\:|\-|\*/, "_")}.vcf.gz.csi") into vcfidx_ch

    script:
    """
    freebayes \
        --region ${interval} \
        --fasta-reference ${fasta} \
        ${params.options} \
        ${aln.collect { "--bam $it" }.join(" ")} \
        | bgzip -c > ${params.project}_${interval.replaceAll(~/\:|\-|\*/, "_")}.vcf.gz
    bcftools index ${params.project}_${interval.replaceAll(~/\:|\-|\*/, "_")}.vcf.gz
    """
}

process make_vcf_list {
    input:
    file(vcf) from makelist_ch.collect()

    output:
    file("vcfs.txt") into vcftxt_ch

    script:
    template 'makelist.py'
}


process merge_vcfs {
    publishDir path: "${params.outdir}/freebayes"

    input:
    file(vcf) from vcf_ch.collect()
    file(vcfidx) from vcfidx_ch.collect()
    file(fasta)
    file(faidx)
    file(vcftxt) from vcftxt_ch

    output:
    file("${params.project}.vcf.gz")
    file("${params.project}.vcf.gz.tbi")

    script:
    """
    bcftools merge -m all --force-samples -l $vcftxt --threads 2 -O z -o ${params.project}_dirty.vcf.gz
    gsort ${params.project}_dirty.vcf.gz $faidx | vcfuniq \
        | bgzip -c > ${params.project}_dirty_sorted.vcf.gz
    bcftools norm -c all -f $fasta --multiallelics - --threads ${task.cpus} \
        --output ${params.project}.vcf.gz --output-type z \
        ${params.project}_dirty_sorted.vcf.gz
    tabix -p vcf ${params.project}.vcf.gz
    """
}
