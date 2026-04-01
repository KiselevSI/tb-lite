#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { TBLITE } from './workflows/tblite'

WorkflowMain.initialise(workflow, params, log)

workflow {
    TBLITE()
}
