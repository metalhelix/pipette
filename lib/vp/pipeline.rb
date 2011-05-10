$:.unshift(File.dirname(__FILE__))

require 'gatk'
require 'annotator'
require 'vcf_filter'

class Pipeline
  def initialize options = {}
    @options = options
    @gatk = GATK.new @options
    # do this now, in case of errors
    @annotator = nil
    if @options[:annotate]
      @annotator = Annotator.new @options
    end
  end

  def execute command
    report command
    #result = %x[#{command}]
    result = system(command)
    #report result
  end

  def report status
    puts "#{Time.now} - " + status if @options[:verbose]
  end

  def realign input_filename, output_prefix
    report "Starting realignment interval creation"
    intervals_file = "#{output_prefix}.intervals"
    params = {"-T" => "RealignerTargetCreator",
              "-I" => input_filename,
              "-o" => intervals_file}

    @gatk.execute params


    report "Starting realignment using intervals"
    realign_bam_file = "#{output_prefix}.realigned.bam"
    params = {"-T" => "IndelRealigner",
              "-I" => input_filename,
              "-targetIntervals" => intervals_file,
              "-o" => realign_bam_file}

    @gatk.execute params

    report "Starting index of realigned BAM with Samtools"
    command = "samtools index #{realign_bam_file}"
    execute command
    realign_bam_file
  end

  def recalibrate input_filename, realign_bam_file, output_prefix
    report "Starting recalibration"
    covar_table_file = "#{output_prefix}.covariate_table.csv"

    params = {"-T" => "CountCovariates",
              "-I" => input_filename,
              "-recalFile" => covar_table_file,
              "-cov" => ["ReadGroupcovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate"]}
    @gatk.execute params

    recal_bam_file = "#{output_prefix}.realigned.recal.bam"

    params = {"-T" => "TableRecalibration",
              "-I" => realign_bam_file,
              "-recalFile" => covar_table_file,
              "--out" => recal_bam_file}
    @gatk.execute params

    cleaned_bam_file = "#{prefix}.realigned.recal.sorted.bam"
    report "Starting sort and index recalibrated BAM"
    command = "samtools sort #{recal_bam_file} #{cleaned_bam_file}"
    execute(command)
    command = "samtools index #{cleaned_bam_file}"
    execute(command)
    cleaned_bam_file
  end

  def call_snps variant_bam_file, output_prefix
    report "Starting SNP calling"
    snps_vcf_file = "#{output_prefix}.snps.vcf"

    params = {"-T" => "UnifiedGenotyper",
              "-I" => variant_bam_file,
              "-o" => snps_vcf_file,
              "-stand_call_conf" => "30.0",
              "-stand_emit_conf" => "30.0"}
    @gatk.execute params
    snps_vcf_file
  end

  def call_indels variant_bam_file, output_prefix
    report "Starting Indel calling"
    indels_vcf_file = "#{output_prefix}.indels.vcf"

    params = {"-T" => "UnifiedGenotyper",
              "-I" => variant_bam_file,
              "-o" => indels_vcf_file,
              "-glm" => "DINDEL",
              "-stand_call_conf" => "30.0",
              "-stand_emit_conf" => "30.0"}
    @gatk.execute params
    indels_vcf_file
  end

  def filter vcf_files
    report "Filtering vcf files"
    output_files = []
    vcf_files.each do |vcf_file|
      output_file = vcf_file.append_filename ".filtered.vcf"
      VCFFilter.filter vcf_file, output_file
      output_files << output_file
    end
    output_files
  end

  def annotate vcf_files
    report "Annotating vcf files"
    output_files = []
    vcf_files.each do |vcf_file|
      output_files << @annotator.annotate(vcf_file)
    end
    output_files
  end
end
