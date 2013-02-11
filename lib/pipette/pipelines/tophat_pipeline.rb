require 'pipette/pipeline'
require 'pipette/tools'

class TophatPipeline < Pipeline
  name "tophat"
  description "Alignment using tophat"
  options do
    opts = OptionParser.new do |o|
      o.on('-n', '--name SAMPLE_NAME', String, 'REQUIRED - Human readable sample name') {|b| options[:name] = b}
      o.on('-i', '--input FASTQ_FILE', String, 'REQUIRED - Input Fastq file to align') {|b| options[:input] = b.split(",")}
      options[:pair] = []
      o.on('-p', '--pair FASTQ_FILE', String, 'paired-end fastq file to align') {|b| options[:pair] = b.split(",")}
      o.on('-r', '--reference FA_FILE', 'REQUIRED - Reference Fasta file for genome') {|b| options[:reference] = b}
      options[:gtf] = ""
      o.on('-g', '--gtf GTF_FILE', 'GTF file to use, if any.') {|b| options[:gtf] = b}

      o.on('-o', '--output PREFIX', 'Output prefix to use for generated files') {|b| options[:output] = b}
      options[:threads] = 1
      o.on('-t', '--threads NUM', Integer, 'Number of threads to use') {|b| options[:threads] = b}

      options[:tophat] = %x[which tophat].chomp
      o.on('--tophat BIN_PATH', String, "Specify location of tophat") {|b| options[:tophat] = b}
      options[:samtools] = %x[which samtools].chomp
      o.on('--samtools BIN_PATH', String, "Specify location of samtools") {|b| options[:samtools] = b}
      o.on('--picard BIN_PATH', String, "Specify location of picard") {|b| options[:picard] = b}
      options[:tophat_params] = ""
      o.on('--tophat_params BIN_PATH', String, "Specify additional tophat params") {|b| options[:tophat_params] = b}

      o.on('-y', '--yaml YAML_FILE', String, 'Yaml configuration file that can be used to load options. Command line options will trump yaml options') {|b| options.merge!(Hash[YAML::load(open(b)).map {|k,v| [k.to_sym, v]}]) }
      o.on('-s', "--steps STEPS" , Array, 'Specify only which steps of the pipeline should be executed') {|b| options[:steps] = b.collect {|step| step} }
      o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
    end
    opts
  end

  step :check_input do
    run do |inputs, outputs|
      raise "ERROR - input FASTQ file required. Use -i parameter, or -h for more info" unless inputs[:input]
      inputs[:input].each do |input|
        raise "ERROR Input file not found at:#{inputs[:input]}" unless File.exists? input
      end

      if inputs[:pair]
        inputs[:pair].each do |pair|
          raise "ERROR - pair FASTQ file required. Use -p parameter, or -h for more info" unless File.exists? pair
        end
      end
      raise "ERROR - reference Fasta file required. Use -r parameter or -h for more info" unless inputs[:reference]
      raise "ERROR Reference file not found at:#{inputs[:reference]}." unless File.exists? inputs[:reference]
      raise "ERROR samtools not found at:#{inputs[:samtools]}." unless File.exists? inputs[:samtools]
      raise "ERROR tophat not found at:#{inputs[:tophat]}." unless inputs[:tophat] and File.exists? inputs[:tophat]
      report "checking inputs complete"
    end
  end

  step :create_output do
    input :output
    run do |inputs, outputs|
      report "creating output directory if needed"
      output_dir = File.expand_path(File.dirname(inputs[:output]))
      command = "mkdir -p #{output_dir}"
      execute command
    end
  end

  step :align do
    input :name
    input :input
    input :pair
    input :reference
    input :output
    input :tophat
    input :tophat_params
    input :gtf
    input :threads
    output :bam_file do |inputs|
      "#{inputs[:output]}.bam"
    end
    run do |inputs, outputs|
      read_one = inputs[:input].join(",")
      read_two = nil
      read_two = inputs[:pair].join(",") unless inputs[:pair].empty?
      output_prefix = inputs[:output]

      tophat_command = "#{inputs[:tophat]}"
      if inputs[:gtf] and !inputs[:gtf].empty?
        tophat_command += " -G #{inputs[:gtf]}"
      end

      if inputs[:tophat_params]
        tophat_command += " #{inputs[:tophat_params]}"
      end

      tophat_command += " -o #{output_prefix}"
      tophat_command += " -p #{inputs[:threads]}"

      # end of tophat_command - add bowtie index
      tophat_command += " #{inputs[:reference]}"
      # and reads
      tophat_command += " #{read_one}"
      if read_two
        tophat_command += " #{read_two}"
      end

      report tophat_command

      report "Starting tophat align"
      # command = "#{inputs[:bwa]} aln -t #{inputs[:threads]} -o 1 -e 1 #{inputs[:reference]} #{inputs[:input]} > #{output_read_one}"
      # execute command

    end
  end

  step :index_bam do
    input :bam_file
    input :samtools
    output :bam_file do |inputs|
      "#{inputs[:output]}.uniq.sort.rmdup.bam"
    end
    run do |inputs, outputs|
      first_bam_name = "#{inputs[:output]}.uniq.bam"
      command = "#{inputs[:samtools]} view -uS -bq 1 #{inputs[:bam_file]} > #{first_bam_name}"
      # execute command
      second_bam_prefix = first_bam_name.split(".")[0..-2].join(".") + ".sort"
      second_bam_out = second_bam_prefix + ".bam"
      command = "#{inputs[:samtools]} sort #{first_bam_name} #{second_bam_prefix}"
      # execute command
    end
  end

end
