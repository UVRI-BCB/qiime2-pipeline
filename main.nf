#!/usr/bin/env nextflow
/*
========================================================================================
               Q I I M E 2   P I P E L I N E
========================================================================================
              qiime2 NEXTFLOW PIPELINE FOR UVRI
 
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    ============================================================
    galaxyuvri-ea/qiime2-pipeline  ~  version ${params.version}
    ============================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run galaxygalaxyuvri-ea-ea/qiime2-pipeline --reads '*_R{1,2}.fastq' --reference 'silva_classifier.qza' --man manifest.csv
    
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. local / galaxyuvri-ea /umic
      --reference                   Path to taxonomic database to be used for annotation
      --manifestFile                Path to manifest file to be used for importing data into qiime2 artifacts
      --denoiser                    Strategy to be used for denoising data dada2 or deblur 
    
    Other arguments:
      --pool                        Should sample pooling be used to aid identification of low-abundance ASVs? Options are pseudo pooling: "pseudo", true: "T", false: "F"
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    
    Trimming arguments (optional):
      --primerFor                   Forward primer
      --primerRev                   Reverse primer
      --truncFor                    Select minimum acceptable length for R1 (--truncFor). Reads will be truncated at truncFor and reads shorter than this are discarded (default 0, no trimming).
      --truncRev                    Select minimum acceptable length for R2 (--truncRev). Reads will be truncated at truncRev and reads shorter than this are discarded (default 0, no trimming).
      --trimFor                     Set length of R1 (--trimFor) that needs to be trimmed (set 0 if no trimming is needed)
      --trimRev                     Set length of R2 (--trimRev) that needs to be trimmed (set 0 if no trimming is needed)
      --maxEE                       After truncation, R1 reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10)), default=2
      --truncQ                      Truncate reads at the first instance of a quality score less than or equal to truncQ; default=2
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Configurable variables
params.name = false
params.project = false
params.email = false
params.plaintext_email = false

// Show help emssage
params.help = false
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

ref=file(params.reference)
man_file=file(params.man)

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into {data; data2}

// Header log info
log.info "============================================================"
log.info " galaxyuvri-ea/qiime2-pipeline  ~  version ${params.version}"
log.info "============================================================"
def summary = [:]
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Reference']      = params.reference
summary['ManifestFile']   = params.man
summary['Denoiser']       = params.denoiser
summary['Forward primer'] = params.primerFor
summary['Reverse primer'] = params.primerRev
summary['truncFor']       = params.truncFor
summary['trimFor']        = params.trimFor
summary['trimRev']        = params.trimRev
summary['truncFor']       = params.truncFor
summary['truncRev']       = params.truncRev
summary['truncQ']         = params.truncQ
summary['pool']           = params.pool
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) {
    summary['E-mail Address'] = params.email
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
* Run Fastqc 
*/

process runfastQC {
    tag { "runFastqc.${pairId}" }
    publishDir "${params.outdir}/qcresults", mode: "copy", overwrite: false

    input:
    set pairId, file(in_fastq) from data

    output:
    file("${pairId}_fastqc/*.zip") into fastqc_files

    """
    mkdir ${pairId}_fastqc
    fastqc --outdir ${pairId}_fastqc \
    ${in_fastq.get(0)} \
    ${in_fastq.get(1)}
    """
}

process runMultiQC{
    publishDir "${params.outdir}/qcresults", mode: "copy", overwrite: false

    input:
    file('*') from fastqc_files.collect()

    output:
    file('multiqc_report.html')

    """
    multiqc .
    """
}

process importData{
  publishDir("${params.outdir}/imported_data")
  
  label 'big_mem'
 
  input:
  file manFile from  man_file 

  output:
  file "demux.qza" into demux_view, demux_trim, demux_dada2

  """
  qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ${manFile} --output-path demux.qza --input-format PairedEndFastqManifestPhred33
  """
}

/*
process viewImported {
  publishDir "${params.outdir}/imported_data", mode: "copy", overwrite: false

  input:
  file in_data from demux_view

  output:
  file 'demux.qzv'

  """
  qiime demux summarize --i-data ${in_data} --o-visualization demux.qzv
  """
}
*/
process adapter_trimming {
  
  publishDir "${params.outdir}/imported_data", mode: "copy", overwrite: false

  input:
  file in_data from demux_trim

  output:
  file 'adapter_trimmed.qza' into trimmed_view, trimmed_denoise

  """
  qiime cutadapt trim-paired --i-demultiplexed-sequences ${in_data} --p-cores ${params.max_cpus}  --p-front-f ${params.primerFor} --p-front-r ${params.primerRev} --o-trimmed-sequences adapter_trimmed.qza
  """

}
 /* 
process viewTrimmed {
  publishDir "${params.outdir}/imported_data", mode: "copy", overwrite: false

  input:
  file in_data from trimmed_view

  output:
  file 'trimmed.qzv' into visual_demux

  """
  qiime demux summarize --i-data ${in_data} --o-visualization 'trimmed.qzv'
  """
}
*/
if(params.denoiser=='dada2'){

  process run_dada2{

  publishDir "${params.outdir}/denoised/dada2", mode: "copy", overwrite: false

  label 'big_mem'
  
  input:
  file clean_reads from trimmed_denoise

  output:
  file "rep_seq.qza" into seqs
  file "table.qza" into table
  file "deno_stats.qza" into denoised
  
  """
  qiime dada2 denoise-paired --i-demultiplexed-seqs ${clean_reads} \
  --p-trim-left-f ${params.trimFor} \
  --p-trim-left-r ${params.trimRev} \
  --p-trunc-len-f ${params.truncFor} \
  --p-trunc-len-r ${params.truncRev} \
  --p-n-threads ${params.max_cpus} \
  --p-trunc-q ${params.truncQ}\
  --p-chimera-method ${params.pool} \
  --o-table table.qza \
  --o-representative-sequences rep_seq.qza \
  --o-denoising-stats deno_stats.qza 
  """
  }

}else if(params.denoiser=="deblur") {
  
  process merge_paired {

    label 'big_mem'

    publishDir "${params.outdir}/denoised/deblur", mode: "copy", overwrite: false
    
    input:
    file clean_reads from trimmed_denoise

    output:
    file "trimmed_joined.qza" into joined

    """
    qiime vsearch join-pairs --i-demultiplexed-seqs ${clean_reads} --o-joined-sequences trimmed_joined.qza
    """
  }

  process run_deblur {

    label 'big_mem'

    publishDir "${params.outdir}/denoised", mode: "copy", overwrite: false

    input:
    file joined_reads from joined

    output:
    file "rep_seq.qza" into seqs
    file "table.qza" into table
    file "deno_stats.qza" into denoised

    """
    qiime deblur denoise-16S  --i-demultiplexed-seqs ${joined_reads} --p-jobs-to-start ${params.max_cpus} --p-sample-stats --p-trim-length 250 --o-table table.qza --o-representative-sequences rep_seq.qza --o-stats deno_stats.qza
    """

    }
}


process alignment {
  
  publishDir "${params.outdir}/alignment", mode: "copy", overwrite: false

  input:
  file seq from seqs

  output:
  file "aligned_seqs.qza" into aligned
  file "masked_aligned_seqs.qza" into masked

  """
  qiime alignment mafft --i-sequences ${seq} \
  --o-alignment aligned_seqs.qza
  qiime alignment mask --i-alignment aligned_seqs.qza \
  --o-masked-alignment masked_aligned_seqs.qza
  """
}

process tree_construction {
  
  publishDir "${params.outdir}/alignment", mode: "copy", overwrite: false

  input:
  file mask_align from masked

  output:
  file "rooted_tree.qza" into phylo_tree

  """
  qiime phylogeny fasttree --i-alignment ${mask_align} \
  --o-tree unrooted_tree.qza
  qiime phylogeny midpoint-root --i-tree unrooted_tree.qza --o-rooted-tree rooted_tree.qza
  """
}

process classification {

  label 'big_mem'
  
  publishDir "${params.outdir}/taxonomy", mode: "copy", overwrite: false

  input:
  file reference from ref
  file seq from seqs

  output:
  file "tax.qza" into taxa

  """
  qiime feature-classifier classify-sklearn --i-classifier ${reference} --i-reads ${seq}  --o-classification tax.qza
  """
}


process export_artifacts {
    publishDir "${params.outdir}/analysis_report", mode: "copy", overwrite: false

    input:
    file abund from table
    file seq from seqs
    file tax from taxa
    file tree from phylo_tree
    file deno from denoised
    
    output:
    file "feature-table.tsv" into abundance
    file "taxonomy.tsv" into taxonomy
    file "tree.nwk" into tree
    file "stats.tsv" into stats
    file "feature-table.biom" into asv_biom
    file "dna-sequences.fasta" into asv_fasta

    """
    qiime tools export --input-path ${abund} --output-path .
    biom convert -i ./feature-table.biom -o ./feature-table.tsv --to-tsv
    qiime tools export --input-path ${seq} --output-path .	
    qiime tools export --input-path ${tax} --output-path .
    qiime tools export --input-path ${tree} --output-path .
    qiime tools export --input-path ${deno} --output-path .
    """
}

process run_picrust2 {

  label 'big_mem'

  publishDir "${params.outdir}", mode: "copy", overwrite: false

  input:
  file asv_table from asv_biom
  file asv_seq from asv_fasta

  output:
  file "*"

  """
  picrust2_pipeline.py -s ${asv_seq} -i ${asv_table} -o picrust2_out 	
  """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {
  
    def subject = "[galaxyuvri-ea/qiime2-pipeline] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[galaxyuvri-ea/qiime2-pipeline] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
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
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container

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
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[galaxyuvri-ea/qiime2-pipeline] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[galaxyuvri-ea/qiime2-pipeline] Sent summary e-mail to $params.email (mail)"
        }
    }
}
