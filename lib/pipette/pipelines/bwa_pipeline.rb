require 'pipette/pipeline'
require 'pipette/tools'

class BwaPipeline < Pipeline
  name "bwa"
  options do
    opts = OptionParser.new do |o|
      o.on('-i', '--input FASTQ_FILE', 'REQUIRED - Input Fastq file to align') {|b| options[:input] = b}
      o.on('-p', '--pair FASTQ_FILE', String, 'REQUIRED - paired-end fastq file to align') {|b| options[:pair] = b}
      o.on('-r', '--reference FA_FILE', 'REQUIRED - Reference Fasta file for genome') {|b| options[:reference] = b}

      o.on('-o', '--output PREFIX', 'Output prefix to use for generated files') {|b| options[:output] = b}
      o.on('-n', '--name NAME', String, 'REQUIRED - Human readable sample name') {|b| options[:name] = b}

      options[:samtools] = %x[which samtools].chomp
      o.on('--samtools BIN_PATH', String, "Specify location of samtools") {|b| options[:samtools] = b}
      options[:bwa] = %x[which bwa].chomp
      o.on('--bwa BIN_PATH', String, "Specify location of bwa") {|b| options[:bwa] = b}
      o.on('--picard BIN_PATH', String, "Specify location of picard") {|b| options[:picard] = b}

      o.on('-y', '--yaml YAML_FILE', String, 'Yaml configuration file that can be used to load options. Command line options will trump yaml options') {|b| options.merge!(Hash[YAML::load(open(b)).map {|k,v| [k.to_sym, v]}]) }
      o.on('-s', "--steps STEPS" , Array, 'Specify only which steps of the pipeline should be executed') {|b| options[:steps] = b.collect {|step| step} }
      o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
    end
    opts
  end

  step :check_input do
    run do |inputs, outputs|
      raise "ERROR - input FASTQ file required. Use -i parameter, or -h for more info" unless inputs[:input]
      raise "ERROR - pair FASTQ file required. Use -p parameter, or -h for more info" unless inputs[:pair]
      raise "ERROR Input file not found at:#{inputs[:input]}" unless File.exists? inputs[:input]
      raise "ERROR - reference Fasta file required. Use -r parameter or -h for more info" unless inputs[:reference]
      raise "ERROR Reference file not found at:#{inputs[:reference]}." unless File.exists? inputs[:reference]
      raise "ERROR samtools not found at:#{inputs[:samtools]}." unless File.exists? inputs[:samtools]
      raise "ERROR bwa not found at:#{inputs[:bwa]}." unless inputs[:bwa] and File.exists? inputs[:samtools]
      raise "ERROR picard not found at:#{inputs[:picard]}." unless inputs[:picard] and File.exists? inputs[:picard]
      report "creating output directory if needed"
      output_dir = File.expand_path(File.dirname(inputs[:output]))
      command = "mkdir -p #{output_dir}"
      execute command
      report "checking inputs complete"
    end
  end

  step :align do
    input :input
    input :pair
    input :reference
    input :output
    input :bwa
    output :sam_file do |inputs|
      "#{inputs[:output]}.aligned.sam.gz"
    end
    run do |inputs, outputs|
      read_one = inputs[:input]
      read_two = inputs[:pair]
      output_prefix = inputs[:output]

      output_read_one = "#{output_prefix}.read1.sai"
      output_read_two = "#{output_prefix}.read2.sai"

      report "Starting bwa align - first read"
      command = "#{inputs[:bwa]} aln -o 1 -e 1 #{inputs[:reference]} #{inputs[:input]} > #{output_read_one}"
      execute command

      report "Staring bwa align - second read"
      command = "#{inputs[:bwa]} aln -o 1 -e 1 #{inputs[:reference]} #{inputs[:pair]} > #{output_read_two}"
      execute command

      #clean_id = inputs[:id].strip.downcase.gsub(" ", "")

      report "Starting bwa sampe"
      #read_group = "@RG\tID:#{clean_id}\tSM:ds\tPL:Illumina"
      command = "#{inputs[:bwa]} sampe #{inputs[:reference]} #{output_read_one} #{output_read_two} #{inputs[:input]} #{inputs[:pair]} | gzip > #{outputs[:sam_file]}"
      execute command
    end
  end

  step :create_bam do
    input :sam_file
    input :samtools
    output :bam_file do |inputs|
      "#{inputs[:output]}.uniq.sort.rmdup.bam"
    end
    run do |inputs, outputs|
      first_bam_name = "#{inputs[:output]}.uniq.bam"
      command = "#{inputs[:samtools]} view -uS -bq 1 #{inputs[:sam_file]} > #{first_bam_name}"
      execute command
      second_bam_prefix = first_bam_name.split(".")[0..-2].join(".") + ".sort"
      second_bam_out = second_bam_prefix + ".bam"
      command = "#{inputs[:samtools]} sort #{first_bam_name} #{second_bam_prefix}"
      execute command
      command = "#{inputs[:samtools]} rmdup #{second_bam_out} #{outputs[:bam_file]}"
      execute command
      command = "#{inputs[:samtools]} index #{outputs[:bam_file]}"
      execute command
    end
  end

  step :group_bam do
    input :bam_file
    input :name
    input :picard
    input :reference
    output :grouped_bam_file do |inputs|
      bam_prefix = inputs[:bam_file].split(".")[0..-2].join(".")
      "#{bam_prefix}.group.bam"
    end
    run do |inputs, outputs|
      clean_name = inputs[:name].strip.downcase.gsub(" ", "")
      command = "java -jar #{inputs[:picard]}/AddOrReplaceReadGroups.jar INPUT=#{inputs[:bam_file]}"
      command += " OUTPUT=#{outputs[:grouped_bam_file]} SORT_ORDER=coordinate RGLB=1 RGPL=illumina RGPU=1 RGSM=#{clean_name} VALIDATION_STRINGENCY=LENIENT"
      execute command
    end
  end

  step :order_bam do
    input :bam_file
    input :name
    input :picard
    input :reference
    input :grouped_bam_file
    output :reordered_bam_file do |inputs|
      bam_prefix = inputs[:bam_file].split(".")[0..-2].join(".")
      "#{bam_prefix}.group.reorder.bam"
    end
    run do |inputs, outputs|
      clean_name = inputs[:name].strip.downcase.gsub(" ", "")
      command = "java -jar #{inputs[:picard]}/ReorderSam.jar INPUT=#{inputs[:grouped_bam_file]}"
      command += " OUTPUT=#{outputs[:reordered_bam_file]} REFERENCE=#{inputs[:reference]} VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"
      execute command
    end
  end
end
