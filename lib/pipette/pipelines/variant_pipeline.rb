require 'pipette/pipeline'
require 'pipette/tools'

class VariantPipeline < Pipeline

  options do
    opts = OptionParser.new do |o|
      o.on('-i', '--input BAM_FILE', 'REQUIRED - Input BAM file to call SNPs on') {|b| options[:input] = b}
      o.on('-r', '--reference FA_FILE', 'REQUIRED - Reference Fasta file for genome') {|b| options[:reference] = b}
      o.on('-o', '--output PREFIX', 'Output prefix to use for generated files') {|b| options[:output] = b}
      o.on('-j', '--cores NUM', "Specify number of cores to run GATK on. Default: #{options[:cores]}") {|b| options[:cores] = b.to_i}
      o.on('-c', '--recalibrate COVARIATE_FILE', "If provided, recalibration will occur using input covariate file. Default: recalibration not performed") {|b| options[:recalibrate] = b}
      o.on('-a', '--annotate GENOME', 'Annotate the SNPs and Indels using Ensembl based on input GENOME. Example Genome: FruitFly') {|b| options[:annotate] = b}
      o.on('--gatk JAR_FILE', String, "Specify GATK installation") {|b| options[:gatk] = b}
      o.on('--snpeff JAR_FILE', String, "Specify snppEff Jar location") {|b| options[:snpeff] = b}
      o.on('--snpeff_config CONFIG_FILE', String, "Specify snppEff config file location") {|b| options[:snpeff_config] = b}
      o.on('--samtools BIN_PATH', String, "Specify location of samtools") {|b| options[:samtools] = b}
      o.on('-q', '--quiet', 'Turn off some output') {|b| options[:verbose] = !b}
      o.on('-s', "--steps STEPS" , Array, 'Specify only which steps of the pipeline should be executed') {|b| options[:steps] = b.collect {|step| step} }
      o.on('-y', '--yaml YAML_FILE', String, 'Yaml configuration file that can be used to load options. Command line options will trump yaml options') {|b| options.merge!(Hash[YAML::load(open(b)).map {|k,v| [k.to_sym, v]}]) }
      o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
    end
    opts
  end

  step :check_input do
    run do |inputs, outputs|
      raise "ERROR - input BAM file required. Use -i parameter, or -h for more info" unless inputs[:input]
      raise "ERROR - reference Fasta file required. Use -r parameter or -h for more info" unless inputs[:reference]
      raise "ERROR GATK JAR not found at:#{inputs[:gatk]}." unless File.exists? inputs[:gatk]
      raise "ERROR samtools not found at:#{inputs[:samtools]}." unless File.exists? inputs[:samtools]
      raise "ERROR Reference file not found at:#{inputs[:reference]}." unless File.exists? inputs[:reference]
      raise "ERROR Input file not found at:#{inputs[:input]}" unless File.exists? inputs[:input]
      if inputs[:recalibrate]
        raise "ERROR covariate file not found at:#{inputs[:recalibrate]}." unless File.exists? inputs[:recalibrate]
      end
    end
  end

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
    input :filtered_vcf_files
    output :annotation_files do |input|
      #TODO: this is the only step that cannot specify
      #what the output file names will be given the input.
      #perhaps use the results of the run and add them to
      #the later inputs as well? then wouldn't need this
      #block on outputs
      outputs = input[:filtered_vcf_files].collect {|vcf| "#{vcf}.annotation.txt"}
      outputs
    end
    run do |inputs, outputs|
      vcf_files = inputs[:filtered_vcf_files]
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

