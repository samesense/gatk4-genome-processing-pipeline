version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data processing according to the GATK Best Practices (June 2016)
## for human whole-genome and exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

#import "./Alignment.wdl" as Alignment
#import "./SplitLargeReadGroup.wdl" as SplitRG
#import "./Qc.wdl" as QC
#import "./BamProcessing.wdl" as Processing
#import "./Utilities.wdl" as Utils
#import "../structs/GermlineStructs.wdl" as Structs

import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/tasks/Alignment.wdl" as Alignment
import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/tasks/SplitLargeReadGroup.wdl" as SplitRG
import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/tasks/Qc.wdl" as QC
import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/tasks/BamProcessing.wdl" as Processing
import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/tasks/Utilities.wdl" as Utils
import "https://raw.githubusercontent.com/samesense/gatk4-genome-processing-pipeline/1.0.0/structs/GermlineStructs.wdl" as Structs

# WORKFLOW DEFINITION
workflow UnmappedBamToAlignedBam {

  input {
    SampleAndUnmappedBams sample_and_unmapped_bams
    GermlineSingleSampleReferences references
    PapiSettings papi_settings

    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu

    String cross_check_fingerprints_by
    File? haplotype_database_file
    Float lod_threshold
    String recalibrated_bam_basename
  }

  Float cutoff_for_large_rg_in_gb = 20.0

  String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"

  Int compression_level = 2

  # Get the version of BWA to include in the PG record in the header of the BAM produced
  # by MergeBamAlignment.
  call Alignment.GetBwaVersion

  # Get the size of the standard reference files as well as the additional reference files needed for BWA

  # Align flowcell-level unmapped input bams in parallel
  scatter (unmapped_bam in sample_and_unmapped_bams.flowcell_unmapped_bams) {

    Float unmapped_bam_size = size(unmapped_bam, "GiB")

    String unmapped_bam_basename = basename(unmapped_bam, sample_and_unmapped_bams.unmapped_bam_suffix)

    # QC the unmapped BAM
    call QC.CollectQualityYieldMetrics as CollectQualityYieldMetrics {
      input:
        input_bam = unmapped_bam,
        metrics_filename = unmapped_bam_basename + ".unmapped.quality_yield_metrics",
        preemptible_tries = papi_settings.preemptible_tries
    }

    if (unmapped_bam_size > cutoff_for_large_rg_in_gb) {
      # Split bam into multiple smaller bams,
      # map reads to reference and recombine into one bam
      call SplitRG.SplitLargeReadGroup as SplitRG {
        input:
          input_bam = unmapped_bam,
          bwa_commandline = bwa_commandline,
          bwa_version = GetBwaVersion.bwa_version,
          output_bam_basename = unmapped_bam_basename + ".aligned.unsorted",
          reference_fasta = references.reference_fasta,
          compression_level = compression_level,
          preemptible_tries = papi_settings.preemptible_tries
      }
    }

    if (unmapped_bam_size <= cutoff_for_large_rg_in_gb) {
      # Map reads to reference
      call Alignment.SamToFastqAndBwaMemAndMba as SamToFastqAndBwaMemAndMba {
        input:
          input_bam = unmapped_bam,
          bwa_commandline = bwa_commandline,
          output_bam_basename = unmapped_bam_basename + ".aligned.unsorted",
          reference_fasta = references.reference_fasta,
          bwa_version = GetBwaVersion.bwa_version,
          compression_level = compression_level,
          preemptible_tries = papi_settings.preemptible_tries
      }
    }

    File output_aligned_bam = select_first([SamToFastqAndBwaMemAndMba.output_bam, SplitRG.aligned_bam])

    Float mapped_bam_size = size(output_aligned_bam, "GiB")

    # QC the aligned but unsorted readgroup BAM
    # no reference as the input here is unsorted, providing a reference would cause an error
    call QC.CollectUnsortedReadgroupBamQualityMetrics as CollectUnsortedReadgroupBamQualityMetrics {
      input:
        input_bam = output_aligned_bam,
        output_bam_prefix = unmapped_bam_basename + ".readgroup",
        preemptible_tries = papi_settings.preemptible_tries
    }
  }

  # Sum the read group bam sizes to approximate the aggregated bam size
  call Utils.SumFloats as SumFloats {
    input:
      sizes = mapped_bam_size,
      preemptible_tries = papi_settings.preemptible_tries
  }

  # MarkDuplicates and SortSam currently take too long for preemptibles if the input data is too large
  Float gb_size_cutoff_for_preemptibles = 110.0
  Boolean data_too_large_for_preemptibles = SumFloats.total_size > gb_size_cutoff_for_preemptibles

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call Processing.MarkDuplicates as MarkDuplicates {
    input:
      input_bams = output_aligned_bam,
      output_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = sample_and_unmapped_bams.base_file_name + ".duplicate_metrics",
      total_input_size = SumFloats.total_size,
      compression_level = compression_level,
      preemptible_tries = if data_too_large_for_preemptibles then 0 else papi_settings.agg_preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = CollectQualityYieldMetrics.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = CollectUnsortedReadgroupBamQualityMetrics.insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = CollectUnsortedReadgroupBamQualityMetrics.insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_metrics

    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    File output_bam = MarkDuplicates.output_bam
  }
}
