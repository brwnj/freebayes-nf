# Why

Yes, `freebayes-parallel` does work, but will only parallelize across a single node.
Ideally, and what we're doing here, is allowing Nextflow to handle what jobs are
split across -- whether it's a node or multiple nodes.

# Usage

You provide `freebayes-nf` with some intervals, like the .fai of your reference and
specify a width. Sub-intervals are created of width across your original intervals
upon which `freebayes` will operate.

The resultants VCFs are merged and the final VCF is decomposed and normalized using
`bcftools`.

In practice this looks like:

```
nextflow run brwnj/freebayes-nf -latest -resume -profile docker \
    --alignments '*.cram' \
    --fasta human.fasta
```

## Required

+ `--alignments`
    + Aligned sequences in .bam and/or .cram format. Indexes (.bai/.crai) must be present.
+ `--fasta`
    + Reference FASTA. Index (.fai) must exist in same directory.

## Options

+ `--outdir`
    + Base results directory for output.
    + Default: '/.results'
+ `--project`
    + File prefix for merged and annotated VCF files.
    + Default: 'variants'
+ `--width`
    + The genomic window size per variant calling job.
    + Default: 5000000
+ `--options`
    + Arguments to be passed to freebayes command in addition to those already supplied like `--bam`, `--region`, and `--fasta-reference`.
    + Single quote these when specifying on the command line, e.g. --options '--pooled-discrete'.
    + Default: '--pooled-continuous --pooled-discrete --genotype-qualities --report-genotype-likelihood-max --allele-balance-priors-off --min-alternate-fraction 0.03 --min-repeat-entropy 1 --min-alternate-count 2'
+ `--intervals`
    + Picard-style intervals file to use rather than intervals defined in .fai.
    + Something like Broad's interval lists work here if you want to omit masked regions.
        + See: https://software.broadinstitute.org/gatk/download/bundle
        + And use the wgs_calling_regions file for your genome build.
