#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/influenzangs
========================================================================================
 nf-core/influenzangs Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/influenzangs
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/influenzangs nextflow run main.nf --readdir "readpath"
    --refdir "referencepath" --extension "_L001_R{1,2}_001.fastq" --outdir "outputpath"
    --linkdir "linkpath" -profile singularity

    Mandatory arguments:
    --readdir [path]              Path to input data
    --refdir	[path]			        Path to reference files
    --extension [str]             Extension of input fastq files
    --outdir	[path]					    Path to which output is saved
    -profile  [str]               Configuration profile to use. Can use multiple (comma separated)
													        Available: conda, docker, singularity, awsbatch, test

    Options:
    --consensus_by_seg [str]      If set, then consensus outputs are organized by segment rather than id

    Other options:
    --email [str]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
    --max_memory [str 128.GB]	    Set the maximum allocated memory
    --max_cpus [int 16]		        Set the maximum allocated cpus
    --max_time [str 48.h]         Set the maximum time allowed for the whole pipeline
    -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

		AWSBatch options:
		--awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
		--awsregion                   The AWS Region for your AWS Batch job to run on
		""".stripIndent()
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

/*
 * Create a channel for input read files
 */

// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']            = params.reads
summary['References']       = params.genome
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-influenzangs-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/influenzangs Workflow Summary'
    section_href: 'https://github.com/nf-core/influenzangs'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publishDirMode,
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo "${workflow.manifest.version}" &> v_pipeline.txt 2>&1 || true
    echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true
    fastqc --version > v_fastqc.txt 2>&1 || true
    trim_galore --version &> v_trim_galore.txt 2>&1 || true
    samtools --version &> v_samtools.txt 2>&1 || true
    bwa &> v_bwa.txt 2>&1 || true
    picard MarkDuplicates --version &> v_picard_mark_dup.txt 2>&1 || true
    picard CreateSequenceDictionary --version &> v_picard_dict.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt 2>&1 || true
    gatk -version > v_gatk.txt 2>&1 || true
    python --version > v_python.txt 2>&1 || true
    varscan > v_varscan.txt 2>&1 || true
    bcftools --version >v_bcftools.txt 2>&1 || true
    bedtools --version >v_bedtools.txt 2>&1 || true
    clustalw --version >v_clustal.txt 2>&1 || true
    R --version > v_R.txt 2>&1 || true
    faidx --version > v_pyfaidx.txt 2>&1 || true
    bbmap.sh --version > v_bbmap.txt 2>&1 || true
    bc -version > v_bc.txt 2>&1 || true
    bam > v_bamutils.txt 2>&1 || true
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 1 - FastQC
 */

Channel.fromFilePairs(params.reads)
     .into { raw_reads_fastqc; raw_reads_trimgalore }

process fastqc {
    tag "$id"

    publishDir "${params.outdir}/fastqc",  mode: params.publishDirMode,
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(id), file(reads) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

/*
 * STEP 2 - Trim_Galore
 */

process trim_samples {
    tag "$id"

    publishDir "${params.outdir}/trimmed_samples",  mode: params.publishDirMode,
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(id), file(reads) from raw_reads_trimgalore

    output:
    set val(id), file("${id}*val_1.fq"), file("${id}*val_2.fq")  into (trimmed_samples_ch,sample_count_ch)
    file "*_fastqc.{zip,html}" into trimqc_results
    file "*_trimming_report.txt" into trim_reports

    script:
    """
    trim_galore --dont_gzip --fastqc --paired ${reads[0]} ${reads[1]}
    """
}

/*
 * STEP 3 - Index Reference
 */
Channel.fromPath(params.genome)
    .map{path -> [path.simpleName,path]}
    .into{ labelled_reference_ch; reference_tosplit_ch }

process index_reference {
    tag "$reference"

    publishDir "${params.outdir}/reference_indices", mode: params.publishDirMode

  	input:
  	set val(reference),file(referencepath) from labelled_reference_ch

  	output:
    set val(reference),file("*.{amb,ann,bwt,pac,sa}") into bwa_ref_ch
    set val(reference),file(referencepath),file("${reference}.{fasta.fai,dict}") into gatk_mutect_reference_ch
  	"""
  	bwa index -p ${reference} -a is ${referencepath}

  	samtools faidx ${referencepath}

  	picard CreateSequenceDictionary R=${referencepath} O=${reference}.dict
    """
}

/*
 * STEP 4 - BWA Alignment
 */

//Combine nets a cartesian product so eg. [bam1,bam2,bam3,..].combine[ref1,ref2,ref3,...]
//is equal to [[bam1,ref1],[bam1,ref2],...,[bam2,ref1],...]

//Final output is [id, trimmed fastq 1, trimmed fastq 2, reference, bwaindices]
bwa_mem_input_ch = trimmed_samples_ch.combine(bwa_ref_ch)

process build_bwa_mem_alignments {
    tag { id + "-" + reference }

    publishDir "${params.outdir}/bwa_alignment", mode: params.publishDirMode

    input:
    set val(id), file(trimreads1), file(trimreads2), val(reference), file(bwaindices) from bwa_mem_input_ch

    output:
    set val(id), val(reference), file("${id}.${reference}.bam.bai"),
    file("${id}.${reference}.bam") into (built_bam_ch, bwa_count_ch)

    script:
    """
    bwa mem -t 4 -M -B 2 -R "@RG\\tID:${id}\\tLB:${id}\\tSM:${id}\\tPL:ILLUMINA" \\
    ${reference} ${trimreads1} ${trimreads2}| samtools sort -@8 -O BAM -o ${id}.${reference}.bam -

    samtools index ${id}.${reference}.bam
    """
}

/*
 * STEP 5 - Picard duplicate removal
 */

process picard_remove_duplicates{
    tag { id + "-" + reference }

    publishDir "${params.outdir}/bam_no_duplicates", mode: params.publishDirMode,
      saveAs: {filename -> filename.indexOf(".metrics") > 0 ? "metrics/$filename" : "$filename"}

    input:
    set val(id),val(reference),file(bamindex),file(bamfile) from built_bam_ch

    output:
    set val(reference), val(id), file("${id}.${reference}.dedup.bam.bai"),
    file("${id}.${reference}.dedup.bam") into bam_no_dup_ch

    file "${id}.${reference}.dedup.metrics" into picardmetrics

    script:
 	  """
 	  picard MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT \\
 		INPUT=$bamfile OUTPUT=${id}.${reference}.dedup.bam \\
 		METRICS_FILE=${id}.${reference}.dedup.metrics TMP_DIR=tmp

    samtools index ${id}.${reference}.dedup.bam
    """
 }


/*
 * STEP 6 - GATK Mutect
 */

bam_no_dup_ch.combine(gatk_mutect_reference_ch, by: 0)
    .into{ GATK_input_mutect_ch; GATK_input_realign_ch }

process GATK_mutect {
    tag {id + "-" + reference}

    publishDir "${params.outdir}/gatk_mutect", mode: params.publishDirMode

    input:
    set val(reference), val(id), file(dupindex), file(dupfile),
    file(reffasta), file(refindex) from GATK_input_mutect_ch

    output:
    set val(id),val(reference),file("${id}.${reference}.prelim.vcf"),
    file("${id}.${reference}.prelim.vcf.idx") into GATK_output_mutect_ch

    script:
    """
    gatk --java-options -Xmx32g Mutect2 -R ${reffasta} -I ${dupfile} --tumor-sample ${id} \\
    -O ${id}.${reference}.prelim.vcf
    """
}

/*
 * STEP 7 - GATK Indel Realignment
 */

process GATK_realign {
    tag {id + "-" + reference}

    label 'OldversionGATK'

    publishDir "${params.outdir}/gatk_realign", mode: params.publishDirMode

    input:
    set val(reference), val(id), file(dupindex), file(dupfile),
    file(reffasta), file(refindex) from GATK_input_realign_ch

    output:
    set val(id), val(reference), file("${id}.${reference}.pe.realigned.bam"),
    file("${id}.${reference}.pe.realigned.bai"), file(reffasta),file(refindex) into GATK_output_realign_ch

    script:
    """
    java -Xmx32g -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reffasta} \\
 		-I ${dupfile} -o ${id}.${reference}.intervals

    java -Xmx32g -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner -R ${reffasta} \\
 		-I ${dupfile} -targetIntervals ${id}.${reference}.intervals \\
 		-o ${id}.${reference}.pe.realigned.bam
    """
 }

/*
 * STEP 8 - GATK Base Recalibration
 */

GATK_input_recalib_ch = GATK_output_mutect_ch.combine(GATK_output_realign_ch, by:[0,1])

process GATK_recalib {
    tag {id + "-" + reference}

    label 'OldversionGATK'

    publishDir "${params.outdir}/GATK_recalib", mode: params.publishDirMode

    input:
    set val(id),val(reference),file(mutectvcf),file(mutectidx),file(realignbam),
    file(realignidx),file(reffasta),file(refindex) from GATK_input_recalib_ch

    output:
    set val(id),val(reference), file("${id}.${reference}.pe.recalib.bam"),
    file("${id}.${reference}.pe.recalib.bai") into (GATK_output_recalib_ch, gatk_count_ch, gatk_extract_fastq_ch)

    script:
    """
    java -Xmx32g -jar /usr/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${reffasta} \\
 		-I ${realignbam} -knownSites ${mutectvcf} \\
 		-o ${id}.${reference}.recalibration.table

    java -Xmx32g -jar /usr/GenomeAnalysisTK.jar -T PrintReads -R ${reffasta} \\
 		-I ${realignbam} -BQSR ${id}.${reference}.recalibration.table \\
 		-o ${id}.${reference}.pe.recalib.bam
    """
 }

/*
 * STEP 9 - Split bam by flu segment
 */

process split_bam {
    tag {id + "-" + reference}

    publishDir "${params.outdir}/alignments_split/${id}/${reference}", mode: params.publishDirMode

    input:
    set val(id),val(reference), file(gatkbam), file(gatkbai) from GATK_output_recalib_ch

    output:
    file("${id}.${reference}_*.bam") optional true into split_bam_ch

    script:
    """
    bam splitChromosome --in ${gatkbam}  --out ${id}. --bamout
    """
}

/*
 * STEP 10 - Split fasta files by flu segment
 */

process split_reference {
    tag "$reference"

    publishDir "${params.outdir}/split_ref", mode: params.publishDirMode

    input:
    set val(reference),file(referencepath) from reference_tosplit_ch

    output:
    file("${reference}_*.fasta") into split_reference_ch

    script:
    """
    faidx -x $referencepath
    """
}

/*
 * STEP 11 - Get fasta sizes
 */

split_reference_ch
    .flatten()
 		.map{path -> [path.simpleName,path]}
    .into{ split_ref_ch; fastasize_input_ch; ref_make_cons_ch; ref_align_1_ch; ref_align_pre_ch; ref_dummy_align_ch }

process split_fasta_size {
    tag "$reference"

    publishDir "${params.outdir}/split_fastasize", mode: params.publishDirMode

    input:
    set val(reference),file(referencepath) from fastasize_input_ch

    output:
 	  set val(reference),file("${reference}.fasta.sizes") into (split_fasta_size_ch, fasta_size_align_ch)

 	  """
 	  picard CreateSequenceDictionary R=${referencepath} O=${reference}.dict

 	  grep "^@SQ" ${reference}.dict | awk '{ print \$2"\t"\$3}' | \\
 	  sed "s;..:;;g" > ${reference}.fasta.sizes
 	  """
 }

/*
 * STEP 12 - Make Coverage Summary
 */

//Recover keys from file name of split bam files. After mapping we get the format [reference_segment, sample id, bam file]
//Combining with split_ref_ch nets us [reference_segment, sample id, bam file, reference segment fasta file]
split_bam_ch
 		.flatten()
 		.map{file -> [file.getFileName().toString().split("\\.")[1], file.getFileName().toString().split("\\.")[0], file]}
 		.combine(split_ref_ch, by:0)
    .into{ mpileup_input_ch; cov_input_nosize_ch; bam_choose_consensus_ch; bcftools_bam_ch }

//Combining with split_fasta_size_ch give us [reference_segment, sample id, bam file, reference segment fasta file, fasta.sizes file]
 cov_input_ch = cov_input_nosize_ch
 		.combine(split_fasta_size_ch, by:0)

process build_coverage {
    tag {id + "-" + refseg}

    publishDir "${params.outdir}/coverage_sum", mode: params.publishDirMode

    input:
    set val(refseg), val(id), file(splitbam), file(refseg_fasta), file(fastasizes) from cov_input_ch

    output:
    set val(id),val(refseg),file("${id}.${refseg}.cov.summary") into (cov_summary_choose_consensus_ch, cov_summary_plot_ch, big_cov_input_ch)
    set val(id),val(refseg),file("${id}.${refseg}.cov") into (cov_choose_consensus_ch, cov_graph_ch)

    script:
    """
    bedtools genomecov -d -ibam ${splitbam} -g  ${fastasizes} | grep ${refseg} > \\
  	${id}.${refseg}.cov

    coverage_summary_from_genomecov.perl -i ${id}.${refseg}.cov -b ${fastasizes} \\
 		-o ${id}.${refseg}.cov.summary
 	  """
 }

/*
 * STEP 13 - Make mpileup
 */

process make_mpileup {
    tag {id + "-" + referenceseg}

    publishDir "${params.outdir}/mpileup", mode: params.publishDirMode

    input:
    set val(referenceseg), val(id), file(splitbam), file(ref_fasta) from mpileup_input_ch

    output:
    set val(id), val(referenceseg), file("${splitbam}.mpileup") into (mpileup_make_consensus_ch, mpileup_varscan_ch)

    script:
    """
    samtools mpileup -a -f ${ref_fasta} -O -d 100000 -B -o \\
    ${splitbam}.mpileup ${splitbam}
    """
}

/*
 * STEP 13.5 - Generate inter-reference alignments
 */

ref_align_2_ch = ref_align_pre_ch
    .join(fasta_size_align_ch, by:0)
    .map{refseg, fastafile, fastasizes-> [refseg.split("_")[1], fastafile, fastasizes]}

reference_to_align_ch = ref_align_1_ch
    .map{refseg, fastafile -> [refseg.split("_")[1], fastafile]}
    .combine(ref_align_2_ch, by:0)
    .filter { it[1] != it[2]}

process align_references {
    tag {shred_fasta + "-" + ref_fasta}

    publishDir "${params.outdir}/reference_cov", mode: params.publishDirMode

    input:
    set val(segment), file(shred_fasta), file(ref_fasta),
    file(ref_fasta_size) from reference_to_align_ch

    output:
    set file("*.*.cov"), val(segment) into aligned_references_true_ch

    script:
    """
    shred_id=\$(echo "${shred_fasta}" | cut -f 1 -d '.')
    ref_id=\$(echo "${ref_fasta}" | cut -f 1 -d '.')

    randomreads.sh -Xmx640m ref=${shred_fasta} out=\${shred_id}.fastq reads=100000 adderrors=f paired=t minlength=300 maxlength=350 illuminanames=t

    bwa index ${ref_fasta}
    bwa mem -p -t 4 -M -B 2 ${ref_fasta} \${shred_id}.fastq | samtools sort -@8 -O BAM -o \${shred_id}.\${ref_id}.bam -

    bedtools genomecov -d -ibam \${shred_id}.\${ref_id}.bam -g ${ref_fasta_size} > \${shred_id}.\${ref_id}.cov
    """
}

process dummy_align_references {
    tag "$refseg"

    publishDir "${params.outdir}/reference_cov", mode: params.publishDirMode

    input:
    set val(refseg), file(fastafile) from ref_dummy_align_ch

    output:
    set file("${refseg}.${refseg}.cov"), val(refseg) into dummy_align_out_ch

    script:
    """
    echo "dummy" > ${refseg}.${refseg}.cov
    """
}

aligned_references_ch = dummy_align_out_ch
  .map{covfile, refseg -> [covfile, refseg.split("_")[1]]}
  .concat(aligned_references_true_ch)

/*
 * STEP 14 - Choose best reference
 */

//Group the individual channel entries by influenza segments across all references
//The map function extracts segment information from the reference segment
//GroupTuple groups all entries by the same id and segments
//Final channel contents (per entry) are as follows:
//[id,segment,list of bam/coverage summary/pileup files, list of reference fastas (only bam channel)]
grouped_bam_choose_ch = bam_choose_consensus_ch
    .map{refseg,id,bam,refsegfasta -> [id,refseg.split("_")[1],bam,refsegfasta]}
    .groupTuple(by:[0,1])

grouped_cov_summary_choose_ch = cov_summary_choose_consensus_ch
    .map{id,refseg,covsum -> [id,refseg.split("_")[1],covsum]}
    .groupTuple(by:[0,1])

grouped_cov_choose_ch = cov_choose_consensus_ch
    .map{id,refseg,covsum -> [id,refseg.split("_")[1],covsum]}
    .groupTuple(by:[0,1])

grouped_aligned_references_ch = aligned_references_ch
      .groupTuple(by: 1)

//Combine the above 3 channels into a single channel with content as follows
//[id, segment, bam file list, reference fasta list, mpileup file list, summary file list]
choose_consensus_input_ch = grouped_bam_choose_ch
    .combine(grouped_cov_summary_choose_ch, by:[0,1])
    .combine(grouped_cov_choose_ch, by:[0,1])
    .combine(grouped_aligned_references_ch, by: 1)

process choose_best_references {
    tag {id + "-" + segment}

    publishDir "${params.outdir}/finalchosenref", mode: params.publishDirMode, pattern: "*.finalchosen"
    publishDir "${params.outdir}/pass_fail_mix", mode: params.publishDirMode, pattern: "*.filter"

    input:
    set val(segment), val(id), file(bamlist), file(fastalist), file(summarylist), file(covlist), file(referencecov) from choose_consensus_input_ch

    output:
    file("${id}.${segment}.finalchosen") into (finalchosen_ch,plotchosen_ch)
    file("${id}.${segment}.filter") into chosenfilter_ch

    script:
    """
    get_best_reference.R ${id} ${segment}
    """
}

/*
 * STEP 15 - Generate consensus sequence
 */

finalchosen_ch
    .splitCsv()
    .flatten()
    .map{text -> [text.split("\\.")[0], text.split("\\.")[1], text.split("\\.")[2]]}
    .into{ consensus_keys_ch ; chosen_drop_filter_ch}

chosen_drop_filter_ch
    .map{ id, refseg, state -> [id, refseg] }
    .into{varscan_keys_ch; bcftools_keys_ch; bamdepth_keys_ch}

//Final structure is {reference_segment, sample id, state, mpileupfile, reference_segment fasta}
make_consensus_input_ch = consensus_keys_ch
    .combine(mpileup_make_consensus_ch, by: [0,1])
    .map{ id,refseg, state, mpileupfile -> [refseg,id,state,mpileupfile] }
    .combine(ref_make_cons_ch, by:0 )

process make_consensus {
    tag {id + "-" + chosenrefseg}

    publishDir "${params.outdir}/consensus/$state", mode: params.publishDirMode,
      saveAs: {
        filename -> filename.indexOf(".log") > 0 ? "log/$filename"
        : params.consensus_by_seg ? "${chosenrefseg.split("_")[1]}/$filename"
        : "${id}/$filename"
        }

    input:
    set val(chosenrefseg), val(id), val(state), file(mpileupfile), file(chosenfasta) from make_consensus_input_ch

    output:
    set val(id), val(state), file("${id}.${chosenrefseg}.fa") optional true into consensus_to_concat_ch
    file "${mpileupfile}.log"

    script:
    """
    build_consensus_from_variants.perl -i ${mpileupfile} -r ${chosenfasta} -l 0.2 \\
    -u 0.8 -c 10 -o ${id}.${chosenrefseg}.fa
    """
}

/*
 * STEP - Concatenate all consensus all_segments
 */

grouped_consensus_to_concat_ch =  consensus_to_concat_ch
    .filter { it[1] == "PASS" }
    .map{ id,state,file -> [id,file] }
    .groupTuple()

process concatenate_consensus {
    tag "$id"

    publishDir "${params.outdir}/consensus", mode: params.publishDirMode

    input:
    set val(id), file(consfiles) from grouped_consensus_to_concat_ch

    output:
    file "${id}.cons.fasta"

    script:
    """
    cat ${consfiles} > ${id}.cons.fasta
    """
}

/*
 * STEP 16 - Generate filter report
 */

process create_filter_sheet {
    publishDir "${params.outdir}", mode: params.publishDirMode

 	  input:
 	  file filter from chosenfilter_ch.collect()

 	  output:
 	  file "all_samples_filtered.csv" into filter_report

 	  """
 	  filter_sheet.R
 	  """
}

/*
 * STEP 17 - Variant calling with varscan
 */
chosen_varscan_ch = varscan_keys_ch.combine(mpileup_varscan_ch, by:[0,1])

process call_varscan {
    tag {id + "-" + refseg}

    publishDir "${params.outdir}/prevariants", mode: params.publishDirMode

    input:
    set val(id), val(refseg), file(mpileupfile) from chosen_varscan_ch

    output:
    set val(id),file("${id}.${refseg}.varscan.vcf.{gz,gz.tbi}") optional true into vcf_varscan_ch

    script:
    """
    varscan mpileup2cns ${mpileupfile} --variants --min-avg-qual 30 --min-var-freq 0.01 \\
    --output-vcf 1 > ${id}.${refseg}.varscan.vcf
    if [ \$(cat ${id}.${refseg}.varscan.vcf | grep -v '^#' | wc -c) -ne 0 ]; then
        bgzip ${id}.${refseg}.varscan.vcf
        tabix ${id}.${refseg}.varscan.vcf.gz
    fi
    """
}

/*
 * STEP 18 - Variant calling with bcftools
 */

chosen_bcftools_ch = bcftools_keys_ch
    .map{id,refseg -> [refseg,id]}
    .combine(bcftools_bam_ch, by:[0,1])

process call_bcftools {
    tag {id + "-" + refseg}

    publishDir "${params.outdir}/prevariants", mode: params.publishDirMode

    input:
    set val(refseg), val(id), file(recalibfile), file(refsegfasta) from chosen_bcftools_ch

    output:
    set val(id), file("${id}.${refseg}.multiallelic.bcftools.vcf.{gz,gz.tbi}") optional true into vcf_bcftools_ch

    script:
    """
    bcftools mpileup -Ou -a AD --per-sample-mF --redo-BAQ -d 10000 --min-BQ 20 -f ${refsegfasta} ${recalibfile} | bcftools call -P 1.0e-2 -m -v -Ov | \\
    vcfutils.pl varFilter -a 1 -w 1 -W 3 -d 5 -1 0.05 -2 0.05 -3 0.05 -4 0.05 > ${id}.${refseg}.multiallelic.bcftools.vcf
    if [ \$(cat ${id}.${refseg}.multiallelic.bcftools.vcf | grep -v '^#' | wc -c) -ne 0 ]; then
        bgzip ${id}.${refseg}.multiallelic.bcftools.vcf
        tabix ${id}.${refseg}.multiallelic.bcftools.vcf.gz
    fi
    """
 }

/*
 * STEP 19 - Concatenate all segments
 */

varscan_combine_ch = vcf_varscan_ch
    .groupTuple()
    .map { id, file -> [id, file.flatten(), "varscan"]}

varscan_bcftools_ch = vcf_bcftools_ch
    .groupTuple()
    .map{ id, file -> [id, file.flatten(), "bcftools"]}

vcf_concat_input_ch = varscan_combine_ch
    .concat(varscan_bcftools_ch)

process vcf_concat {
    tag {id + "-" + vcf_label}

    publishDir "${params.outdir}/prevariants", mode: params.publishDirMode

    input:
    set val(id), file(vcf_files), val(vcf_label) from vcf_concat_input_ch

    output:
    set val(id), file("${id}.${vcf_label}.concat.vcf.gz"), file("${id}.${vcf_label}.concat.vcf.gz.tbi") into vcf_concat_ch

    script:
    """
    vcf-concat ${id}.*.vcf.gz > ${id}.${vcf_label}.concat.vcf
    bgzip ${id}.${vcf_label}.concat.vcf
    tabix ${id}.${vcf_label}.concat.vcf.gz
    """
}

/*
 * STEP 20 - Combine varscan and bcftools variant files
 */

process combine_vcf{
    tag "$id"

    publishDir "${params.outdir}/variants", mode: params.publishDirMode

  	input:
  	set val(id), file(vcf_file), file(vcf_index) from vcf_concat_ch.groupTuple()

  	output:
  	set val(id), file("${id}.variants.vcf") into combine_vcf_ch

  	script:
  	"""
  	vcf-merge ${vcf_file} > ${id}.variants.vcf
  	"""
}

/*
 * STEP 21 - Get Variant frequency tables
 */

 process get_variant_frequency {
     tag "$id"

     publishDir "${params.outdir}/variant_frequencies", mode: params.publishDirMode

     input:
     set val(id), file(vcf_file) from combine_vcf_ch

     output:
     file "${id}.variants.csv"

     script:
     """
     short_vcf.R ${vcf_file} > ${id}.variants.csv
     """
 }

/*
 * STEP 22 - Plot 10x coverage
 */

groupedplotchosen_ch = plotchosen_ch
    .map{file->[file.getFileName().toString().split("\\.")[1],file]}
    .groupTuple()

Rplotin_ch = cov_summary_plot_ch
    .map{id,refseg,covsum -> [refseg.split("_")[1],covsum]}
    .groupTuple()
    .combine(groupedplotchosen_ch, by: 0)

process plot_graphs {
    tag "$segment"

    publishDir "${params.outdir}/reports/coverage_graphs", mode: params.publishDirMode

    input:
    set val(segment), file(coverage_sum), file(chosenref) from Rplotin_ch

    output:
    file("10x_coverage_${segment}.png")
    file("10x_coverage_${segment}.html")

    script:
    """
    plotscript_perseg.R ${segment} cov.summary
    """
}

/*
 * STEP 23 - Get big coverage table
 */

process get_big_table {
    tag "$id"

    publishDir "${params.outdir}/reports/coverage_tables", mode: params.publishDirMode

    input:
    set val(id), val(refseg), file(coverage_sum) from big_cov_input_ch.groupTuple()

    output:
    file("${id}.all_coverage.csv")

    script:
    """
    create_big_coverage.R ${id} cov.summary
    """
}

/*
 * STEP 24 - Get readcounts from bam alignments
 */

count_input_ch = bwa_count_ch
		.combine(gatk_count_ch,by:[0,1])
		.groupTuple( by:0,sort:true )
		.combine(sample_count_ch, by:0)

process count_report {
    tag "$id"

    publishDir "${params.outdir}/reports/count_reports", mode: params.publishDirMode

    input:
	  set val(id), val(referencelist), file(bwaindex), file(bwafile),
    file(gatkfile), file(gatkindex), file(read1),file(read2) from count_input_ch

	  output:
	  file("counts.${id}.xls") into count_table_ch

    script:
    bashreflist = referencelist.collect{"$it"}.join(' ')

    """
    echo -ne "\t${id}\n" > counts.${id}.xls
    COUNT_READ1=`echo \$( cat ${read1}|wc -l )/4|bc`
    COUNT_READ2=`echo \$( cat ${read2}|wc -l )/4|bc`
    COUNTS=\$(( COUNT_READ1 + COUNT_READ2 ))
    echo -ne "\t\${COUNTS}\n" >> counts.${id}.xls

    for reference in ${bashreflist}
    do
        echo -ne "\${reference} RAW\t" >> counts.${id}.xls

		    BAM_COUNTS=`samtools view -F 260 -c \\
        ${id}.\${reference}.bam`

		    echo -ne "\${BAM_COUNTS}\n" >> counts.${id}.xls

		    echo -ne "\${reference} RECAL\t" >> counts.${id}.xls

		    RECALIB_COUNTS=`samtools view -F 260 -c \\
		    ${id}.\${reference}.pe.recalib.bam`

	      echo -ne "\${RECALIB_COUNTS}\n" >> counts.${id}.xls
    done
	"""
}

/*
 * STEP 25 - Combine all count tables
 */

process make_list {
    publishDir "${params.outdir}/reports", mode: params.publishDirMode

    input:
    file counttablelist from count_table_ch.collect(sort:true)

    output:
    file "all_samples.counts.xls"

    script:
    """
    COMMAND=(paste)
    for counttable in ${counttablelist}
    do
		    COMMAND+=( \${counttable} )
    done
    "\${COMMAND[@]}" > all_samples.counts.xls
	"""
}

/*
 * STEP 26 - Extract Flu reads-only fastq files
 */

process extract_fastq {
    tag "$id"

    publishDir "${params.outdir}/extract_fastq", mode: params.publishDirMode

    input:
    set val(id), val(reference), file(gatklist), file(gatkindex) from gatk_extract_fastq_ch.groupTuple()

	  output:
	  file("${id}.R1.fastq.gz")
    file("${id}.R2.fastq.gz")

    script:
    """
    for gatkfile in ${gatklist}
    do
        samtools view -b -F 4 \${gatkfile} > \${gatkfile}.filtered
    done

    samtools merge ${id}.merged.bam *.bam.filtered

    samtools sort -n ${id}.merged.bam | samtools fastq -N -F 256 -1 \\
		${id}.R1.fastq.gz -2 \\
		${id}.R2.fastq.gz -
	"""
}

depth_input_ch = bamdepth_keys_ch
    .combine(cov_graph_ch, by:[0,1])

/*
 * STEP 28 - Combine depth files by sample id and plot
 */

process plot_depth{
    tag "$id"

    publishDir "${params.outdir}/depth_plots", mode: params.publishDirMode

    input:
    set val(id), val(refseg), file(depthfiles) from depth_input_ch.groupTuple(sort:true)

    output:
    file("${id}.basedepth.pdf")

    script:
    """
    for depthfile in $depthfiles
    do
        cat \${depthfile} >> ${id}.all_segments.depth
    done

    pandocloc=`which pandoc`

    base_cov.R ${id}.all_segments.depth ${id} \${pandocloc}
    """
}



/*
 * STEP 27 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: params.publishDirMode

    input:
    file multiqc_config from ch_multiqc_config
    file mqc_custom_config from ch_multiqc_custom_config.collect().ifEmpty([])
    // TODO nf-core: Add in log files from your new processes for MultiQC to find!
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('trimqc/*') from trimqc_results.collect().ifEmpty([])
    file ('trimrep/*') from trim_reports.collect().ifEmpty([])
    file ('picardmetrics/*') from picardmetrics.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''

    """
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}

/*
 * STEP 28 - Output Description HTML
 */

process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/influenzangs] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/influenzangs] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nf-core/influenzangs] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/influenzangs] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/influenzangs] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[nf-core/influenzangs] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file

    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/influenzangs]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/influenzangs]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/influenzangs v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
