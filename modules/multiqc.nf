/*
    MULTIQC
    Input : collected list of all QC output files
    Output: multiqc_report.html and multiqc_data/
*/

process MULTIQC {
    label      'process_single'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --title    "nf-gatk-pipeline QC Report"
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    """
}
