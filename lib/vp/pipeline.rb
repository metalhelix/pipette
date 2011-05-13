$:.unshift(File.dirname(__FILE__))

require 'gatk'
require 'annotator'
require 'vcf_filter'

class Pipeline
  @@valid_steps = [:realign, :recalibrate, :call, :filter, :annotate]

  attr_accessor :prefix, :steps

  def initialize options = {}
    @options = options
    @gatk = GATK.new @options
    # do this now, in case of errors
    @annotator = nil
    if @options[:annotate]
      @annotator = Annotator.new @options
    end
    self.prefix = nil
  end

  def self.valid_steps
    @@valid_steps
  end

  def check_steps
    if !steps || steps.empty?
      raise "ERROR no steps provided"
    end
    steps.each do |step|
      unless @@valid_steps.include? step
        raise "ERROR: #{step} not a valid step.\nvalid steps: #{@@valid_steps.join(",")}"
      end
    end
  end

  def performing_step? step
    @steps.include? step.to_sym
  end

  def execute command
    report command
    result = system(command)
  end

  def report status
    puts "#{Time.now} - " + status if @options[:verbose]
  end

  def realign input_filename, output_prefix
    # output file from realign step
    realign_bam_file = "#{output_prefix}.realigned.bam"

    if performing_step? :realign
      report "Starting realignment interval creation"
      intervals_file = "#{output_prefix}.intervals"
      params = {"-T" => "RealignerTargetCreator",
                "-I" => input_filename,
                "-o" => intervals_file}

      @gatk.execute params

      report "Starting realignment using intervals"

      params = {"-T" => "IndelRealigner",
                "-I" => input_filename,
                "-targetIntervals" => intervals_file,
                "-o" => realign_bam_file}

      @gatk.execute params

      report "Starting index of realigned BAM with Samtools"
      command = "samtools index #{realign_bam_file}"
      execute command
    else
      report "Skipping realignment"
    end
    report "Outputfile: #{realign_bam_file}"
    realign_bam_file
  end

  def recalibrate input_filename, realign_bam_file, output_prefix
    cleaned_bam_file = "#{prefix}.realigned.recal.sorted.bam"
    if performing_step? :recalibrate
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


      report "Starting sort and index recalibrated BAM"
      command = "samtools sort #{recal_bam_file} #{cleaned_bam_file}"
      execute(command)
      command = "samtools index #{cleaned_bam_file}"
      execute(command)
    else
      report "Skipping Recalibration step"
    end

    report "Recalibration output: #{cleaned_bam_file}"
    cleaned_bam_file
  end

  def call_snps variant_bam_file, output_prefix
    snps_vcf_file = "#{output_prefix}.snps.vcf"
    if performing_step? :call
      report "Starting SNP calling"

      params = {"-T" => "UnifiedGenotyper",
                "-I" => variant_bam_file,
                "-o" => snps_vcf_file,
                "-stand_call_conf" => "30.0",
                "-stand_emit_conf" => "30.0"}
      @gatk.execute params
    else
      report "Skipping SNP calling"
    end
    report "SNP Calling output: #{snps_vcf_file}"
    snps_vcf_file
  end

  def call_indels variant_bam_file, output_prefix
    indels_vcf_file = "#{output_prefix}.indels.vcf"
    if performing_step? :call
      report "Starting Indel calling"

      params = {"-T" => "UnifiedGenotyper",
                "-I" => variant_bam_file,
                "-o" => indels_vcf_file,
                "-glm" => "DINDEL",
                "-stand_call_conf" => "30.0",
                "-stand_emit_conf" => "30.0"}
      @gatk.execute params
    else
      report "Skipping Indel calling"
    end
    report "Indel calling output: #{indels_vcf_file}"
    indels_vcf_file
  end

  def filter vcf_files
    output_files = []
    vcf_files.each do |vcf_file|
      output_file = vcf_file.append_filename ".filtered.vcf"
      if performing_step? :filter
        report "Filtering vcf file: #{vcf_file}"
        VCFFilter.filter vcf_file, output_file
      else
        report "Skipping vcf filter of #{vcf_file}"
      end
      output_files << output_file
    end
    report "VCF filter output: #{output_files.join(", ")}"
    output_files
  end

  def annotate vcf_files
    output_files = []
    vcf_files.each do |vcf_file|
      if performing_step? :annotate
        report "Annotating #{vcf_file}"
        output_file = @annotator.annotate(vcf_file)
      else
        report "Skipping annotation of #{vcf_file}"
      end
      output_files << output_file 
    end
    report "annotation output #{output_files.join(", ")}"
    output_files
  end
end

