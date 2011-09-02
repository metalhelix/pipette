require 'pipette/pipeline'
require 'pipette/tools'

class BwaPipeline < Pipeline

  options do
    opts = OptionParser.new do |o|
      o.on('-i', '--input FASTQ_FILE', 'REQUIRED - Input Fastq file to align') {|b| options[:input] = b}
      o.on('-o', '--output PREFIX', 'Output prefix to use for generated files') {|b| options[:output] = b}
      o.on('-r', '--reference FA_FILE', 'REQUIRED - Reference Fasta file for genome') {|b| options[:reference] = b}
      o.on('-p', '--pair FASTQ_FILE', String, 'If present, paired-end fastq file') {|b| options[:pair] = b}
      o.on('-n', '--name NAME', String, 'REQUIRED - Human Readable Name') {|b| options[:name] = b}
      o.on('-g', '--group NAME', String, 'REQUIRED - Group Name') {|b| options[:group] = b}

      o.on('-y', '--yaml YAML_FILE', String, 'Yaml configuration file that can be used to load options. Command line options will trump yaml options') {|b| options.merge!(Hash[YAML::load(open(b)).map {|k,v| [k.to_sym, v]}]) }
      o.on('-s', "--steps STEPS" , Array, 'Specify only which steps of the pipeline should be executed') {|b| options[:steps] = b.collect {|step| step} }
      o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
    end
    opts
  end

  step :align_first do
    input :input
    input :output
    output :align_out do |inputs|
      ["#{inputs[:output]}.sai"]
    end
    run do |inputs, outputs|
      input_filename = inputs[:input]
    end
  end

  step :align_second do
    input :output
    input :align_out
    output :align_out do |inputs|
      inputs[:align_out] << "#{inputs[:output]}_2.sai"
    end
    run do |inputs, outputs|
      puts "align 2"
    end
  end

  step :join do
    run do |inputs, outputs|
    end
  end

  step :create_bam do
    run do |inputs, outputs|
    end
  end

  step :sort do
    run do |inputs, outputs|
    end
  end

  step :remove_dups do
    run do |inputs, outputs|
    end
  end

  step :index do
    run do |inputs, outputs|
    end
  end
end
