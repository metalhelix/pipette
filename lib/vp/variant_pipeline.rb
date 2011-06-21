
require 'vp/pipeline'
require 'vp/gatk'
require 'vp/annotator'
require 'vp/vcf_filter'

class VariantPipeline < Pipeline

  def realign options
    input_filename = options[:input]
    output_prefix = options[:output]
    # output file from realign step
    realign_bam_file = "#{output_prefix}.realigned.bam"

    report "Starting realignment interval creation"
    intervals_file = "#{output_prefix}.intervals"
    params = {"-T" => "RealignerTargetCreator",
              "-I" => input_filename,
              "-o" => intervals_file}

    gatk = GATK.new options
    gatk.execute params

    report "Starting realignment using intervals"

    params = {"-T" => "IndelRealigner",
              "-I" => input_filename,
              "-targetIntervals" => intervals_file,
              "-o" => realign_bam_file}

    gatk.execute params

    report "Starting index of realigned BAM with Samtools"
    samtools_path = options[:samtools]
    raise "ERROR: no samtools at: #{samtools_path}" unless File.exists? samtools_path
    command = "#{samtools_path} index #{realign_bam_file}"
    execute command
    report "Outputfile: #{realign_bam_file}"
    results = {:bam_file => realign_bam_file}
    results
  end

  def recalibrate options
    input_filename = options[:input]
    realign_bam_file = options[:bam_file]
    output_prefix = options[:output]

    cleaned_bam_file = "#{prefix}.realigned.recal.sorted.bam"
    report "Starting recalibration"
    covar_table_file = "#{output_prefix}.covariate_table.csv"

    params = {"-T" => "CountCovariates",
              "-I" => input_filename,
              "-recalFile" => covar_table_file,
              "-cov" => ["ReadGroupcovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate"]}
    gatk = GATK.new options
    gatk.execute params

    recal_bam_file = "#{output_prefix}.realigned.recal.bam"

    params = {"-T" => "TableRecalibration",
              "-I" => realign_bam_file,
              "-recalFile" => covar_table_file,
              "--out" => recal_bam_file}
    gatk.execute params


    report "Starting sort and index recalibrated BAM"
    command = "samtools sort #{recal_bam_file} #{cleaned_bam_file}"
    execute(command)
    command = "samtools index #{cleaned_bam_file}"
    execute(command)

    report "Recalibration output: #{cleaned_bam_file}"
    results = {:bam_file => cleaned_bam_file}
  end

  def call_snps options
    variant_bam_file = options[:bam_file]
    output_prefix = options[:output]

    snps_vcf_file = "#{output_prefix}.snps.vcf"
    report "Starting SNP calling"

    params = {"-T" => "UnifiedGenotyper",
              "-I" => variant_bam_file,
              "-o" => snps_vcf_file,
              "-stand_call_conf" => "30.0",
              "-stand_emit_conf" => "30.0"}
    gatk = GATK.new options
    gatk.execute params
    report "SNP Calling output: #{snps_vcf_file}"
    results = {}
    if options[:vcf_files]
      results[:vcf_files] = options[:vcf_files]
      results[:vcf_files] << snps_vcf_file
    else
      results[:vcf_files] = [snps_vcf_file]
    end
    results
  end

  def call_indels options
    variant_bam_file = options[:bam_file]
    output_prefix = options[:output]
    indels_vcf_file = "#{output_prefix}.indels.vcf"
    report "Starting Indel calling"

    params = {"-T" => "UnifiedGenotyper",
              "-I" => variant_bam_file,
              "-o" => indels_vcf_file,
              "-glm" => "DINDEL",
              "-stand_call_conf" => "30.0",
              "-stand_emit_conf" => "30.0"}
    gatk = GATK.new options
    gatk.execute params
    report "Indel calling output: #{indels_vcf_file}"
    results = {}
    if options[:vcf_files]
      results[:vcf_files] = options[:vcf_files]
      results[:vcf_files] << indels_vcf_file
    else
      results[:vcf_files] = [indels_vcf_file]
    end
    results
  end

  def filter options
    vcf_files = options[:vcf_files]
    output_files = []
    vcf_files.each do |vcf_file|
      output_file = vcf_file.append_filename ".filtered.vcf"
      report "Filtering vcf file: #{vcf_file}"
      VCFFilter.filter vcf_file, output_file
      output_files << output_file
    end
    report "VCF filter output: #{output_files.join(", ")}"
    results = {:vcf_files => output_files}
  end

  def annotate options
    vcf_files = options[:vcf_files]
    output_files = []
    vcf_files.each do |vcf_file|
      report "Annotating #{vcf_file}"
      annotator = Annotator.new
      output_file = annotator.run(vcf_file, options)
      output_files << output_file
    end
    report "annotation output #{output_files.join(", ")}"
    results = {:annotation_files => output_files}
    results
  end
end

