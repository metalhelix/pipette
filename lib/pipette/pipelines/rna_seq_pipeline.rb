require 'pipette/pipeline'
require 'pipette/tools'

class RnaSeqPipeline < Pipeline
  name "rna_seq"
  description "Alignment using tophat. analysis with cufflinks and other tools"
  options do
    opts = OptionParser.new do |o|
      o.on('-n', '--name SAMPLE_NAME', String, 'REQUIRED - Human readable sample name') {|b| options[:name] = b}
      o.on('-i', '--input FASTQ_FILE', Array, 'REQUIRED - Input Fastq file to align') {|b| options[:input] = b}
      options[:pair] = []
      o.on('-p', '--pair FASTQ_FILE', Array, 'paired-end fastq file to align') {|b| options[:pair] = b}
      o.on('--index INDEX_PREFIX', 'REQUIRED - index prefix for bowtie to use') {|b| options[:index] = b}
      o.on('-r', '--reference FA_FILE', 'Reference Fasta file for genome') {|b| options[:reference] = b}
      options[:gtf] = ""
      o.on('-g', '--gtf GTF_FILE', 'GTF file to use, if any.') {|b| options[:gtf] = b}

      o.on('-o', '--output PREFIX', 'Output prefix to use for generated folders') {|b| options[:output] = b}
      options[:threads] = 1
      o.on('-t', '--threads NUM', Integer, 'Number of threads to use') {|b| options[:threads] = b}

      options[:tophat] = %x[which tophat].chomp
      o.on('--tophat BIN_PATH', String, "Specify location of tophat") {|b| options[:tophat] = b}
      options[:cufflinks] = %x[which cufflinks].chomp
      o.on('--cufflinks BIN_PATH', String, "Specify location of cufflinks") {|b| options[:cufflinks] = b}
      options[:samtools] = %x[which samtools].chomp
      o.on('--samtools BIN_PATH', String, "Specify location of samtools") {|b| options[:samtools] = b}
      o.on('--picard BIN_PATH', String, "Specify location of picard") {|b| options[:picard] = b}
      options[:tophat_params] = ""
      o.on('--tophat_params BIN_PATH', String, "Specify additional tophat params") {|b| options[:tophat_params] = b}
      options[:cufflinks_params] = ""
      o.on('--cufflinks_params BIN_PATH', String, "Specify additional cufflinks params") {|b| options[:cufflinks_params] = b}

      o.on('-y', '--yaml YAML_FILE', String, 'Yaml configuration file that can be used to load options. Command line options will trump yaml options') {|b| options.merge!(Hash[YAML::load(open(b)).map {|k,v| [k.to_sym, v]}]) }
      o.on('-s', "--steps STEPS" , Array, 'Specify only which steps of the pipeline should be executed') {|b| options[:steps] = b.collect {|step| step} }
      o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
    end
    opts
  end

  step :check_input do
    run do |inputs, outputs|
      # begin
        raise "ERROR - input FASTQ file required. Use -i parameter, or -h for more info" unless inputs[:input]
        [inputs[:input]].flatten.each do |input|
          raise "ERROR Input file not found at:#{input}" unless File.exists? File.expand_path(input)
        end

        if inputs[:pair]
          inputs[:pair].each do |pair|
            raise "ERROR - pair FASTQ file required. Use -p parameter, or -h for more info" unless File.exists? pair
          end
        end
        raise "ERROR - bowtie index file required. Use --index parameter or -h for more info" unless inputs[:index]
        raise "ERROR samtools not found at:#{inputs[:samtools]}." unless File.exists? File.expand_path(inputs[:samtools])
        raise "ERROR tophat not found at:#{inputs[:tophat]}." unless inputs[:tophat] and File.exists? File.expand_path(inputs[:tophat])
        raise "ERROR cufflinks not found at:#{inputs[:cufflinks]}." unless inputs[:cufflinks] and File.exists? File.expand_path(inputs[:cufflinks])
        report "checking inputs complete"
      # rescue Exception => e
      #   puts e
      #   puts "Error found in check input. Please fix and run again"
      #   exit(1)
      # end
    end
  end


  step :create_tophat_out do
    input :output
    input :name
    output :tophat_output_path do |inputs|
      File.expand_path(File.join(inputs[:output],inputs[:name],"tophat"))
    end
    run do |inputs, outputs|
      report "creating output directory if needed"
      command = "mkdir -p #{outputs[:tophat_output_path]}"
      execute command
    end
  end

  step :tophat do
    input :input
    input :pair
    input :index
    input :tophat
    input :tophat_params
    input :gtf
    input :threads
    input :tophat_output_path
    output :bam_file do |inputs|
      File.join(inputs[:tophat_output_path],"accepted_hits.bam")
    end
    run do |inputs, outputs|
      read_one = [inputs[:input]].flatten.join(",")
      read_two = nil
      read_two = [inputs[:pair]].flatten.join(",") unless inputs[:pair].empty?
      output_prefix = inputs[:tophat_output_path]

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
      tophat_command += " #{inputs[:index]}"
      # and reads
      tophat_command += " #{read_one}"
      if read_two
        tophat_command += " #{read_two}"
      end

      report tophat_command

      report "Starting tophat align"
      execute(tophat_command)

    end
  end

  step :create_cufflinks_out do
    input :output
    input :name
    output :cufflinks_output_path do |inputs|
      File.expand_path(File.join(inputs[:output], inputs[:name], "cufflinks"))
    end
    run do |inputs, outputs|
      report "creating cufflinks directory if needed"
      command = "mkdir -p #{outputs[:cufflinks_output_path]}"
      execute command

    end
  end

  step :cufflinks do
    input :bam_file
    input :cufflinks_output_path
    input :cufflinks
    input :cufflinks_params
    input :threads
    input :gtf
    output :fpkm_file do |inputs|
      File.join(inputs[:cufflinks_output_path], "genes.fpkm_tracking")
    end
    run do |inputs, outputs|
      cufflinks_command = "#{inputs[:cufflinks]}"
      if inputs[:gtf] and !inputs[:gtf].empty?
        cufflinks_command += " -G #{inputs[:gtf]}"
      end

      if inputs[:cufflinks_params]
        cufflinks_command += " #{inputs[:cufflinks_params]}"
      end

      cufflinks_command += " -o #{inputs[:cufflinks_output_path]}"
      cufflinks_command += " -p #{inputs[:threads]}"

      cufflinks_command += " #{inputs[:bam_file]}"

      report "starting cufflinks run"
      execute(cufflinks_command)
    end
  end

  step :create_output_dir do
    input :output
    output :output_bam_dir do |inputs|
      File.join(inputs[:output], "bams")
    end
    output :output_fpkm_dir do |inputs|
      File.join(inputs[:output], "cufflink_fpkms")
    end
    run do |inputs, outputs|
      report "creating output directory if needed"
      command = "mkdir -p #{outputs[:output_bam_dir]}"
      execute command
      command = "mkdir -p #{outputs[:output_fpkm_dir]}"
      execute command
    end
  end

  step :move_bam_file do
    input :bam_file
    input :output_bam_dir
    input :name
    output :moved_bam_file do |inputs|
      File.join(inputs[:output_bam_dir], "#{inputs[:name]}.bam")
    end
    run do |inputs, outputs|
      command = "mv #{inputs[:bam_file]} #{outputs[:moved_bam_file]}"
      execute command
    end
  end

  step :move_fpkm_file do
    input :fpkm_file
    input :output_fpkm_dir
    input :name
    output :moved_fpkm_file do |inputs|
      File.join(inputs[:output_fpkm_dir], "#{inputs[:name]}.fpkm_tracking")
    end
    run do |inputs, outputs|
      command = "mv #{inputs[:fpkm_file]} #{outputs[:moved_fpkm_file]}"
      execute command
    end
  end

  step :index_bam do
    input :moved_bam_file
    input :samtools
    output :bam_file do |inputs|
      "#{inputs[:moved_bam_file]}.bai"
    end
    run do |inputs, outputs|
      output_dir = File.dirname(inputs[:moved_bam_file])
      Dir.chdir(output_dir) do
        command = "samtools index #{File.basename(inputs[:moved_bam_file])}"
        execute command
      end
    end
  end

end
