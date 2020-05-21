#!/usr/bin/env nextflow
/*
========================================================================================
                         czbiohub/sc2-ngs-analysis
========================================================================================
 #### Homepage / Documentation
 https://github.com/czbiohub/sc2-ngs-analysis
----------------------------------------------------------------------------------------
*/

def helpMessage() {
  log.info nfcoreHeader()
  log.info"""
    Pipeline for running augur commands to build auspice visualization from main output.

    Usage:

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
      --sample_sequences            Glob pattern of sequences
      --sample_metadata             TSV of metadata for nextstrain
      --ref                         Reference FASTA file (default: ${params.ref})
      --ref_gb                      Reference Genbank file (default: ${params.ref_gb})
      --blast_sequences             FASTA of sequences for BLAST alignment
      --nextstrain_sequences        FASTA of sequences to build a tree with
      --nextstrain_metadata         TSV from GISAID
      --minLength                   Minimum base pair length to allow assemblies to pass QC (default: ${params.minLength})
      --maxNs                       Max number of Ns to allow assemblies to pass QC (default: ${params.maxNs})

    Optional arguments:
      --clades                      TSV with clades from nextstrain (default: ${params.clades})
      --sample_vcfs                 Glob pattern of corresponding VCF files
      --gisaid_names                TSV to rename samples to public identifiers, need columns sample_name, gisaid_name
      --qpcr_primers                qPCR primer BED file (default: ${params.qpcr_primers})


    Other options:
      --outdir                      The output directory where the results will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      -resume                       Use cached results

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Reference genomes

ref_fasta = file(params.ref, checkIfExists: true)
ref_gb = file(params.ref_gb, checkIfExists: true)

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}

if ( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

/*
 * Create a channel for input assemblies
 */

 Channel
  .fromPath(params.sample_sequences)
  .ifEmpty { exit 1, "Cannot find any files matching: ${params.sample_sequences}\nNB: Path needs to be enclosed in quotes!" }
  .map {file -> tuple(file.simpleName, file) }
  .into {blastconsensus_in; realign_fa; stats_fa}

Channel
  .fromPath(params.sample_sequences)
  .ifEmpty { exit 1, "Cannot find any files matching: ${params.sample_sequences}\nNB: Path needs to be enclosed in quotes!" }
  .set {merge_fastas_ch}


// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Sample Sequences']            = params.sample_sequences
summary['Fasta Ref']        = params.ref
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile == 'awsbatch') {
  summary['AWS Region']     = params.awsregion
  summary['AWS Queue']      = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
  summary['E-mail Address']    = params.email
  summary['E-mail on failure'] = params.email_on_fail
  summary['MultiQC maxsize']   = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'sc2-ngs-analysis-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'czbiohub/sc2-ngs-analysis Workflow Summary'
    section_href: 'https://github.com/czbiohub/sc2-ngs-analysis'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*
 * STEP 1 - Realign consensus genome
 */
process realignConsensus {
    tag { sampleName }
    publishDir "${params.outdir}/realigned-seqs", mode: 'copy'

    input:
    tuple(sampleName, path(in_fa)) from realign_fa
    path(ref_fasta)

    output:
    tuple(sampleName, path("${sampleName}.realigned.bam")) into realigned_bam
    path("${sampleName}.realigned.bam.bai")

    script:
    """
    minimap2 -ax asm5 -R '@RG\\tID:${sampleName}\\tSM:${sampleName}' \
      ${ref_fasta} ${in_fa} |
      samtools sort -O bam -o ${sampleName}.realigned.bam
    samtools index ${sampleName}.realigned.bam
    """
}

// Redirect queue channel into additional channels for callVariants and combinedVariants
realigned_bam.into { call_variants_bam; combined_variants_bams }
combined_variants_bams = combined_variants_bams.map { it[1] }.collect()

/* 
 * STEP 2 - call variants against the reference
 */

// Allow VCF files to be optional in case we want to include assemblies that aren't produced by the core
// consensus calling pipeline
if (params.sample_vcfs) {
  Channel
    .fromPath(params.sample_vcfs)
    .map {file -> tuple(file.simpleName, file) }
    .set {variants_ch}
} 
else {
  process callVariants {
    tag { sampleName }
    publishDir "${params.outdir}/sample-variants", mode: 'copy'

    input:
    tuple(sampleName, path(in_bam)) from call_variants_bam
    path(ref_fasta)

    output:
    tuple(sampleName, path("${sampleName}.vcf")) into variants_ch
    path("${sampleName}.bcftools_stats") into bcftools_stats_ch

    script:
    """
    bcftools mpileup -f ${ref_fasta} ${in_bam} |
      bcftools call --ploidy 1 -m -P ${params.bcftoolsCallTheta} -v - \
      > ${sampleName}.vcf
    bcftools stats ${sampleName}.vcf > ${sampleName}.bcftools_stats
    """
  }
}

// redirect VCFs to assignClades and primer variants
variants_ch.into {primer_variants_ch; assignclades_in; variants_ch}

clades = file(params.clades, checkIfExists: true)

/*
 * STEP 3 - assign clades
 */
process assignClades {
    // Use Nextstrain definitions to assign clades based on mutations

    input:
    tuple(sampleName, path(vcf)) from assignclades_in
    path(ref_gb)
    path(clades)

    output:
    tuple(sampleName, path("${sampleName}.clades")) into assignclades_out

    script:
    """
    assignclades.py \
        --reference ${ref_gb} --clades ${clades} \
        --vcf ${vcf} --sample ${sampleName}
    """
}

// Search qPCR primer regions if provided
if (params.qpcr_primers) {
    qpcr_primers = file(params.qpcr_primers, checkIfExists: true)
} else {
    qpcr_primers = Channel.empty()
}

/*
 * STEP 4 - search qPCR primers
 */
process searchPrimers {
    publishDir "${params.outdir}/primer-variants", mode: 'copy'

    input:
    tuple(sampleName, path(vcf)) from primer_variants_ch
    path(qpcr_primers)

    output:
    tuple(sampleName, path("${sampleName}_primers.vcf")) into primer_variants_vcf
    path("${sampleName}_primers.primer_variants_stats") into primer_stats_ch

    when:
    params.qpcr_primers

    script:
    """
    bgzip ${vcf}
    bcftools index ${vcf}.gz
    bcftools view -R ${qpcr_primers} -o ${sampleName}_primers.vcf ${vcf}.gz
    bcftools stats ${sampleName}_primers.vcf > ${sampleName}_primers.primer_variants_stats
    """
}

/*
 * STEP 5 - BUILD BLAST DATABASE
 */
// build a BLAST database only if sequences are provided
blast_sequences = params.blast_sequences ? file(params.blast_sequences, checkIfExists: true) : Channel.empty()

process buildBLASTDB {

    input:
    path(blast_sequences)

    output:
    path("blast_seqs.nt*") into blastdb_ch
    path("blast_sequences.fasta") into blast_clean_ch

    script:
    """
    normalize_gisaid_fasta.sh ${blast_sequences} blast_seqs_clean.fasta
    sed 's/\\.//g' blast_seqs_clean.fasta > blast_sequences.fasta

    makeblastdb -in blast_sequences.fasta -parse_seqids -title 'blastseqs' -dbtype nucl -out blast_seqs.nt
    """
}

/*
 * STEP 6 - BLAST CONSENSUS AND GET CLOSEST SEQUENCES
 */
process blastConsensus {
    tag {sampleName}
    publishDir "${params.outdir}/samples/${sampleName}", mode: 'copy'

    input:
    tuple(sampleName, path(assembly)) from blastconsensus_in
    path(dbsequences) from blast_clean_ch
    path(blastdb) from blastdb_ch
    path(ref_fasta)

    output:
    path("${sampleName}.blast.tsv")
    tuple(sampleName, path("${sampleName}_nearest_blast.fasta")) into (nearest_neighbor, collectnearest_in)

    script:
    """
    get_top_hit.py --minLength ${params.minLength} \
      --sequences ${dbsequences} \
      --sampleName ${sampleName} \
      --assembly ${assembly} \
      --default ${ref_fasta}
    """
}

/*
 * STEP 7 - COLLECT ALL THE SEQUENCES
 */
process collectNearest {
    publishDir "${params.outdir}/BLAST", mode: 'copy'

    input:
    path(fastas) from collectnearest_in.map{it[1]}.collect()

    output:
    path("included_samples.fasta") into (included_fastas_ch, contextual_fastas_ch)

    script:
    """
    cat ${fastas} > all_included_samples.fasta
    seqkit rmdup all_included_samples.fasta > deduped_included_samples.fasta
    seqkit faidx -f -r deduped_included_samples.fasta '^[^MN908947.3].*' > included_samples.fasta
    """
}

/*
 * STEP 8 - COMPUTE PER SAMPLE STATISTICS
 */
// combine all the relevant stats files
stats_fa
  .join(variants_ch)
  .join(primer_variants_vcf)
  .join(assignclades_out)
  .join(nearest_neighbor)
  .set{stats_ch_in}

process computeStats {
    tag { sampleName }
    publishDir "${params.outdir}/coverage-plots", mode: 'copy',
        saveAs: { x -> x.endsWith(".png") ? x : null }

    input:
    tuple(sampleName,
          file(in_fa),
          file(vcf),
          file(primer_vcf),
          file(in_clades),
          file(neighbor_fasta)) from stats_ch_in

    output:
    path("${sampleName}.stats.json") into stats_ch

    script:
    """
    alignment_assembly_stats.py \
        --sample_name ${sampleName} \
        --assembly ${in_fa} \
        --vcf ${vcf} \
        --primervcf ${primer_vcf} \
        --neighborfasta ${neighbor_fasta} \
        --clades ${in_clades} \
        --out_prefix ${sampleName}
    """
}

/*
 * STEP 9 - COMBINE VARIANT FILES
 */
process combinedVariants {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(in_bams) from combined_variants_bams
    path(ref_fasta)

    output:
    path("combined.vcf") into combined_variants_vcf

    script:
    """
    bcftools mpileup -f ${ref_fasta} ${in_bams} |
      bcftools call --ploidy 1 -m -P ${params.bcftoolsCallTheta} -v - \
      > combined.vcf
    """
}

/*
 * STEP 10 - MERGE ALL ASSEMBLIES
 */
 process mergeAllAssemblies {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(in_fasta) from merge_fastas_ch.collect()

    output:
    path("combined.fa") into merged_assemblies_ch

    script:
    """
    cat ${in_fasta} > combined.fa
    """
}

/*
 * STEP 11 - MERGE ASSEMBLY STATS
 */
process mergeAssemblyStats {
    publishDir "${params.outdir}/run_analysis-stats", mode: 'copy'

    input:
    path(in_json) from stats_ch.collect()

    output:
    path("combined.stats.tsv") into merged_stats_ch

    script:
    """
    merge_stats.py analysis ${in_json} > combined.stats.tsv
    """
}

/*
 * STEP 12 - FILTER ASSEMBLIES
 */
process filterAssemblies {
    publishDir "${params.outdir}", mode: 'copy',
      saveAs: {x -> x.endsWith(".tsv") ? "run_analysis-stats/$x" : x}

    input:
    path(merged_stats) from merged_stats_ch
    path(merged_assemblies) from merged_assemblies_ch
    path(vcf) from combined_variants_vcf

    output:
    path("filtered.stats.tsv")
    path("filtered.fa") into nextstrain_ch
    path("filtered.vcf")

    script:
    """
    filter_assemblies.py \
        --vcf ${vcf} \
        --max_n ${params.maxNs} --min_len ${params.minLength} \
        --stats ${merged_stats} --fasta ${merged_assemblies} \
        --out_prefix filtered
    """
}

/*
 * STEP 13 - RENAME TO GISAID NAMES IF PROVIDED
 */
gisaid_names = params.gisaid_names ? file(params.gisaid_names, checkIfExists: true) : Channel.empty()
process renameAssemblies {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(in_fa) from nextstrain_ch
  path(gisaid_names)

  output:
  path("renamed_sequences.fasta") into renamed_asm

  when:
  params.gisaid_names

  script:
  """
  rename_samples.py --in_fa ${in_fa} --gisaid_names ${gisaid_names}
  """
}
if (!params.gisaid_names) {
  nextstrain_ch.set{renamed_asm}
}

// Put sample and contextual sequences into one channel
renamed_asm
  .mix(contextual_fastas_ch)
  .collect()
  .set {sample_and_contextual_ch}

/*
 * STEP 14 - MAKE NEXTSTRAIN SEQUENCE FILE
 */
nextstrain_metadata = file(params.nextstrain_metadata, checkIfExists: true)
nextstrain_sequences = file(params.nextstrain_sequences, checkIfExists: true)
sample_metadata = file(params.sample_metadata, checkIfExists: true)
include_file = file(params.include_file, checkIfExists: true)

process combineSequences {
    publishDir "${params.outdir}/nextstrain/data", mode: 'copy'

    input:
    path(sample_sequences) from renamed_asm
    path(nextstrain_sequences)
    
    output:
    path("sequences.fasta")

    script:
    //Normalize the GISAID names using Nextstrain's bash script
    """
    cat ${sample_sequences} ${nextstrain_sequences} > all_sequences.fasta
    seqkit rmdup all_sequences.fasta > deduped_sequences.fasta
    normalize_gisaid_fasta.sh deduped_sequences.fasta sequences.fasta ${params.minLength}
    """
}

/*
 * STEP 15 - MAKE NEXTSTRAIN METADATA FILE
 */
process combineMetadata {
    publishDir "${params.outdir}/nextstrain/data", mode: 'copy'

    input:
    path(nextstrain_metadata)
    path(sample_metadata)

    output:
    path("metadata.tsv")

    script:
    """
    combine_metadata.py --nextstrain_metadata ${nextstrain_metadata} \
        --sample_metadata ${sample_metadata}
    """
}

/*
 * STEP 16 - GET CA SEQUENCES
 */
process getCASequences {
  publishDir "${params.outdir}/nextstrain/config", mode: 'copy'

  input:
  path(nextstrain_metadata)
  path(nextstrain_sequences)

  output:
  path("CA_sequences.txt") into CA_sequences

  script:
  String exclude_where = "date='2020' date='2020-01-XX' date='2020-02-XX' date='2020-03-XX' date='2020-04-XX' date='2020-01' date='2020-02' date='2020-03' date='2020-04'"
  """
  augur filter \
  --sequences ${nextstrain_sequences} \
  --metadata ${nextstrain_metadata} \
  --exclude-where ${exclude_where} \
  --exclude-where 'division!=California' \
  --min-length ${params.minLength} \
  --output CA_sequences.fasta

  cat CA_sequences.fasta | grep '>' | awk -F '>' '{print \$2}' > CA_sequences.txt
  """
}

/*
 * STEP 17 - MAKE NEXTSTRAIN INCLUDE FILE
 */
process combineInclude {
    publishDir "${params.outdir}/nextstrain/config", mode: 'copy'

    input:
    path(include_file)
    path(sample_sequences) from renamed_asm
    path(blast_fastas) from included_fastas_ch
    path(CA_sequences)

    output:
    path("include.txt")

    script:
    """
    cat ${blast_fastas} | grep '>' | awk -F '>' '{print \$2}' > blast_fastas.txt
    cat ${sample_sequences} | grep '>' | awk -F '>' '{print \$2}' > sample_sequences.txt
    cat ${include_file} ${CA_sequences} blast_fastas.txt sample_sequences.txt > include.txt
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[czbiohub/sc2-ngs-analysis] Successful: $workflow.runName"
    if (!workflow.success) {
      subject = "[czbiohub/sc2-ngs-analysis] FAILED: $workflow.runName"
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
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp


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
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
          if ( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[czbiohub/sc2-ngs-analysis] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, email_address ].execute() << email_txt
          log.info "[czbiohub/sc2-ngs-analysis] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if (!output_d.exists()) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[czbiohub/sc2-ngs-analysis]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[czbiohub/sc2-ngs-analysis]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  czbiohub/sc2-ngs-analysis v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname(){
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
