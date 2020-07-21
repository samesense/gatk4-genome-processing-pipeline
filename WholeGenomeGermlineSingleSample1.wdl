version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "./tasks/UnmappedBamToAlignedBam1.wdl" as ToBam1
import "./tasks/AggregatedBamQC.wdl" as AggregatedQC
import "./tasks/Qc.wdl" as QC
import "./tasks/BamToCram.wdl" as ToCram
import "./tasks/VariantCalling.wdl" as ToGvcf
import "./structs/GermlineStructs.wdl"

# import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/tasks/UnmappedBamToAlignedBam1.wdl" as ToBam1
# import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/tasks/AggregatedBamQC.wdl" as AggregatedQC
# import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/tasks/Qc.wdl" as QC
# import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/tasks/BamToCram.wdl" as ToCram
# import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/tasks/VariantCalling.wdl" as ToGvcf
# import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/structs/GermlineStructs.wdl"

# WORKFLOW DEFINITION
workflow WholeGenomeGermlineSingleSample1 {

  String pipeline_version = "1.3"

  input {
    SampleAndUnmappedBams sample_and_unmapped_bams
    GermlineSingleSampleReferences references
    PapiSettings papi_settings
    File wgs_coverage_interval_list

    File? haplotype_database_file
    Boolean provide_bam_output = false
    Boolean use_gatk3_haplotype_caller = false
  }

  # Not overridable:
  Int read_length = 250
  Float lod_threshold = -20.0
  String cross_check_fingerprints_by = "READGROUP"
  String recalibrated_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.duplicates_marked.recalibrated"

  call ToBam1.UnmappedBamToAlignedBam {
    input:
      sample_and_unmapped_bams    = sample_and_unmapped_bams,
      references                  = references,
      papi_settings               = papi_settings,

      contamination_sites_ud = references.contamination_sites_ud,
      contamination_sites_bed = references.contamination_sites_bed,
      contamination_sites_mu = references.contamination_sites_mu,

      cross_check_fingerprints_by = cross_check_fingerprints_by,
      haplotype_database_file     = haplotype_database_file,
      lod_threshold               = lod_threshold,
      recalibrated_bam_basename   = recalibrated_bam_basename
  }

  # Outputs that will be retained when execution is complete
  output {
    File duplicate_metrics = UnmappedBamToAlignedBam.duplicate_metrics
    File output_bam = UnmappedBamToAlignedBam.output_bam
  }
}
