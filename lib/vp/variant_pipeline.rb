
require 'vp/pipeline'
require 'vp/gatk'
require 'vp/annotator'
require 'vp/vcf_filter'

class VariantPipeline < Pipeline

  step :realign do
    input :input
    input :output
    input :samtools
    output :bam_file do |inputs|
      "#{inputs[:output]}.realigned.bam"
    end
    run do |inputs, outputs|
      input_filename = inputs[:input]
      output_prefix = inputs[:output]
      # output file from realign step
      realign_bam_file = "#{output_prefix}.realigned.bam"

      report "Starting realignment interval creation"
      intervals_file = "#{output_prefix}.intervals"
      params = {"-T" => "RealignerTargetCreator",
                "-I" => input_filename,
                "-o" => intervals_file}

      gatk = GATK.new inputs
      gatk.execute params

      report "Starting realignment using intervals"

      params = {"-T" => "IndelRealigner",
                "-I" => input_filename,
                "-targetIntervals" => intervals_file,
                "-o" => realign_bam_file}

      gatk.execute params

      report "Starting index of realigned BAM with Samtools"
      samtools_path = inputs[:samtools]
      raise "ERROR: no samtools at: #{samtools_path}" unless File.exists? samtools_path
      command = "#{samtools_path} index #{realign_bam_file}"
      execute command
    end
  end

  step :recalibrate do
    input :input
    input :bam_file
    input :output
    output :bam_file do |input|
      if input[:recalibrate]
        "#{input[:output]}.realigned.recal.sorted.bam"
      else
        input[:bam_file]
      end
    end
    run do |inputs, outputs|
      input_filename = inputs[:input]
      realign_bam_file = inputs[:bam_file]
      output_prefix = inputs[:output]

      cleaned_bam_file = outputs[:bam_file]
      report "Starting recalibration"
      covar_table_file = "#{output_prefix}.covariate_table.csv"

      params = {"-T" => "CountCovariates",
                "-I" => input_filename,
                "-recalFile" => covar_table_file,
                "-cov" => ["ReadGroupcovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate"]}
      gatk = GATK.new inputs
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
    end
  end

  step :call_snps do
    input :bam_file
    input :output
    output :vcf_files do |input|
      vcf_files = input[:vcf_files] ? input[:vcf_files] : []
      vcf_files << "#{input[:output]}.snps.vcf"
      vcf_files
    end
    run do |inputs, outputs|
      variant_bam_file = inputs[:bam_file]
      output_prefix = inputs[:output]

      snps_vcf_file = outputs[:vcf_files][-1]
      report "Starting SNP calling"

      params = {"-T" => "UnifiedGenotyper",
                "-I" => variant_bam_file,
                "-o" => snps_vcf_file,
                "-stand_call_conf" => "30.0",
                "-stand_emit_conf" => "30.0"}
      gatk = GATK.new inputs
      gatk.execute params
      report "SNP Calling output: #{snps_vcf_file}"
    end
  end

  step :call_indels do
    input :bam_file
    input :output
    output :vcf_files do |input|
      vcf_files = input[:vcf_files] ? input[:vcf_files] : []
      vcf_files << "#{input[:output]}.indels.vcf"
      vcf_files
    end
    run do |inputs, outputs|
      variant_bam_file = inputs[:bam_file]
      output_prefix = inputs[:output]
      indels_vcf_file = outputs[:vcf_files][-1]
      report "Starting Indel calling"

      params = {"-T" => "UnifiedGenotyper",
                "-I" => variant_bam_file,
                "-o" => indels_vcf_file,
                "-glm" => "INDEL",
                "-stand_call_conf" => "30.0",
                "-stand_emit_conf" => "30.0"}
      gatk = GATK.new inputs
      gatk.execute params
      report "Indel calling output: #{indels_vcf_file}"
    end
  end

  step :filter do
    input :vcf_files
    output :filtered_vcf_files do |input|
      outputs = input[:vcf_files].collect {|vcf_file| vcf_file.append_filename ".filtered.vcf"}
      outputs
    end

    run do |inputs, outputs|
      vcf_files = inputs[:vcf_files]
      outputs[:filtered_vcf_files].each_with_index do |output_file,index|
        vcf_file = vcf_files[index]
        VCFFilter.filter vcf_file, output_file
      end
    end
  end

  step :annotate do
    input :vcf_files
    output :annotation_files do |input|
      #TODO: this is the only step that cannot specify
      #what the output file names will be given the input.
      #perhaps use the results of the run and add them to
      #the later inputs as well? then wouldn't need this
      #block on outputs
      outputs = input[:vcf_files].collect {|vcf| "#{vcf}.annotation.txt"}
      outputs
    end
    run do |inputs, outputs|
      vcf_files = inputs[:vcf_files]
      output_files = []
      vcf_files.each do |vcf_file|
        report "Annotating #{vcf_file}"
        annotator = Annotator.new
        output_file = annotator.run(vcf_file, inputs)
        output_files << output_file
      end
      report "annotation output #{output_files.join(", ")}"
      results = {:annotation_files => output_files}
      results
    end
  end
end

