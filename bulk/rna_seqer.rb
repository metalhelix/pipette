#!/usr/bin/env ruby

require 'yaml'
require 'optparse'
require 'parallel'

CUR_DIR = File.expand_path(File.dirname(__FILE__))
PIPELINE = File.expand_path(File.join(CUR_DIR, "..", "rna_seq.rb"))

def assert
  raise "assertion failed !" unless yield
end

class Runner
    HELP = <<-HELP
Usage: rna_seqer [TYPE] [ORDER/STARTING_DIR] <OPTIONS>

TYPE can be one of the following:

  order    - Create samples.yml for an order.
             Expects input to be text file with
             flowcell_id[TAB]flowcell_dir
             for all flowcells in the order.

  flowcell - Create samples.yml for flowcell.
             expects input to be flowcell_dir
             for flowcell.

  check    - Checks the validity of a 
             samples.yml file. Input is
             samples.yml file name.
  

  run      - Run pipeline on a particular 
             samples.yml file. Input is 
             samples.yml file name.

  help     - Prints this text, as well as pipeline
             help and exits

  HELP

  # ---
  # constructor
  # handle initialization of
  # class-wide variables
  # ---
  def initialize(options = {})
    @options = options
    @options['sample_report'] ||= 'Sample_Report.csv'
    @options['processes'] ||= 1

    puts @options.inspect
  end

  # ---
  # parse order file
  # order file expected to be
  # FLOWCELL_ID[TAB]FLOWCELL_PATH
  # ---
  def load_order(order_file)
    if !File.exists?(order_file)
      puts "ERROR: order file not found"
      puts "#{order_file}"
      exit(1)
    end
    orders = File.open(order_file, 'r').read.split("\n").collect {|l| l.split("\t")}
  end

  # ---
  # parse sample report file
  # only returns data for samples that have data in the specified directory
  # ---
  def load_sample_report(root_dir, fc_id)
    sample_report_filename = File.join(root_dir, @options['sample_report'])
    if !File.exists?(sample_report_filename)
      puts "ERROR: #{@options['sample_report']} not found in directory"
      puts "#{root_dir}"
      exit(1)
    end

    samples = []
    samples_data = File.open(sample_report_filename, 'r').read.split("\n")
    header = samples_data.shift.split(",").collect {|d| d.strip}
    samples_data.each do |line|
      data = line.split(",").collect {|d| d.strip}
      sample_hash = Hash[header.zip(data)]
      full_path = File.join(root_dir, sample_hash['output'])
      if File.exists? full_path
        sample_hash['root'] = root_dir
        sample_hash['full'] = full_path
        sample_hash['flowcell'] = fc_id
        samples << sample_hash
      end
    end
    samples
  end

  def condense_samples(samples)
    retain = ['full','illumina index', 'read', 'lane' 'flowcell']
    sorted_samples = []
    samples.each do |sample|
      sample.keep_if {|k,v| retain.include?(k)}
      sorted_samples << Hash[sample.sort_by {|k, v| retain.find_index(k)}]
    end
    sorted_samples
  end

  # ---
  # combines data from multiple sample reports
  # ---
  def join_sample_reports(sample_reports)
    keys = @options['merge_samples'] ? ['sample name'] : ['flowcell', 'lane', 'illumina index']
    sample_data = Hash.new {|h,k| h[k] = []}
    sample_reports.each do |sample_report|
      sample_report.each do |sample|
        key = keys.collect {|k| sample[k]}.join("_")
        sample_data[key] << sample
      end
    end
    samples = []
    sample_data.each do |key, data|
      sample = {}
      sample['name'] = key
      first_reads = data.select {|d| d['read'] == '1'}
      first_reads = first_reads.sort {|a,b| "#{a['flowcell']}_#{a['lane']}" <=> "#{b['flowcell']}_#{b['lane']}"}
      second_reads = data.select {|d| d['read'] == '2'}
      second_reads = second_reads.sort {|a,b| "#{a['flowcell']}_#{a['lane']}" <=> "#{b['flowcell']}_#{b['lane']}"}
      sample['input'] = condense_samples(first_reads)
      sample['pair'] = condense_samples(second_reads)
      samples << sample
    end
    samples
  end

  # ---
  # basic checking before running
  # ---
  def check_samples(samples)
    names = samples.collect {|s| s['name']}
    unique_names = names.uniq
    assert { names == unique_names }

    samples.each do |sample|
      if !sample['pair'].empty?
        assert { sample['pair'].size == sample['input'].size }
        pair_lanes = sample['pair'].collect {|s| s['lane']}
        input_lanes = sample['input'].collect {|s| s['lane']}
        assert { pair_lanes == input_lanes }
        pair_flowcells = sample['pair'].collect {|s| s['flowcell']}
        input_flowcells = sample['input'].collect {|s| s['flowcell']}
        assert { pair_flowcells == input_flowcells }
        pair_flowcells = sample['pair'].collect {|s| s['root']}
        input_flowcells = sample['input'].collect {|s| s['root']}
        assert { pair_flowcells == input_flowcells }
      end
    end
  end

  # ---
  # ---
  def load_sample_reports(orders)
    sample_reports = []
    orders.each do |order|
      sample_reports << load_sample_report(order[1], order[0])
    end
    samples = join_sample_reports(sample_reports)
    samples
  end

  # ---
  # ---
  def execute(command)
    puts command
    system(command)
    puts " ---- "
  end

  # ---
  # ---
  def run_pipeline_on_sample(data)
    command = "#{PIPELINE}"
    command += " --name #{data['name']}"
    command += " --input #{data['input']}"
    if data['pair']
      command += " --pair #{data['pair']}"
    end
    command += " #{data['args'].join(" ")}"
    execute(command)
  end

  # ---
  # ---
  def run_pipeline(samples, args)
    processes = @options['processes']
    Parallel.each(samples, :in_processes => processes) do |sample|
      input_string = sample['input'].collect{ |i| i['full']}.join(",")
      pair_string = nil
      if sample['pair'] and !sample['pair'].empty?
        pair_string = sample['pair'].collect{ |i| i['full']}.join(",")
      end

      sample_data = {'name' => sample['name'],
                     'input' => input_string,
                     'pair' => pair_string,
                     'args' => args}

      run_pipeline_on_sample(sample_data)
    end
  end

  # ---
  # ---
  def output_samples_file(samples, args)
    samples.each do |key, value|
      if value.class == String
        value = value.force_encoding 'UTF-8'
      end
    end
    output = {"samples" => samples, "args" => args}
    YAML::ENGINE.yamler='syck'
    File.open(@options['samples_file'],'w') do |file|
      file.write(YAML.dump(output))
    end
  end

  # ---
  # ---
  def parse_samples_file(filename)
    if !File.exists?(filename)
      puts "ERROR: it doesn't look like you gave me"
      puts " a samples.yml file."
      puts "Input file does not exist."
      puts "Input file: #{filename}."
      exit(1)
    end
    config = YAML.load(File.open(filename, 'r'))
    samples = {}
    if config['samples']
      samples = config['samples']
    end
    args = []
    if config['args']
      args = config.delete('args')
    end
    [samples, args]
  end

  # ---
  # ---
  def prepare_order(order_file, args)
    orders = load_order(order_file)
    samples = load_sample_reports(orders)
    output_samples_file(samples, args)
  end

  # ---
  # ---
  def prepare_flowcell(starting_dir, args)
    flowcell_id = File.basename(starting_dir)
    sample_report = load_sample_report(starting_dir, flowcell_id)
    samples = join_sample_reports([sample_report])
    output_samples_file(samples, args)
  end

  def run(input_file, additional_args)

    samples, args = parse_samples_file(input_file)
    args.concat(additional_args)
    run_pipeline(samples, args)
  end

  # ---
  # ---
  def run_sample(args)
    command = "#{PIPELINE} #{args.join(" ")}"
    execute command
  end

  def check(filename)
    samples, args = parse_samples_file(filename)
    check_samples(samples)
  end

  def self.help()
    puts Runner::HELP
    puts system("#{PIPELINE} -h")
  end

  # ---
  # ---
  def self.execute(command, options, args)
    case command
    when "order"
      self.new(options).prepare_order(args.shift, args)
    when "flowcell"
      self.new(options).prepare_flowcell(args.shift, args)
    when "run"
      self.new(options).run(args.shift, args)
    when "check"
      self.new(options).check(args.shift)
    when "sample"
      self.new(options).run_sample(args)
    else
      self.help()
    end
  end
end

options = {}

options['sample_report'] = "Sample_Report.csv"
options['samples_file'] = "./samples.yml"
options['processes'] = 1
options['merge_samples'] = true

opts = OptionParser.new do |o|
  o.on('-y', '--yaml YAML_FILE', String, 'Yaml configuration file that can be used to load options. Command line options will trump yaml options') {|b| options.merge!(Hash[YAML::load(open(b)).map {|k,v| [k.to_sym, v]}]) }
  o.on('--sample_report [SAMPLE_REPORT_NAME]', String, "Name of Sample_Report file. default: #{options['sample_report']}") {|b|  options['sample_report'] = b}
  o.on('--samples_file [SAMPLES_FILE_NAME]', String, "Name of samples output file. default: #{options['samples_file']}") {|b|  options['samples_file'] = b}
  o.on('--processes [NUM_PROCESSES]', Integer, "Number of samples to run at a time. default: #{options['processes']}") {|b|  options['processes'] = b}
  o.on('--no_merge_samples', "Disables merging samples with the same name. default: #{options['merge_samples']}") {|b|  options['merge_samples'] = false}
  o.on('-h', '--help', 'Displays help screen, then exits') {puts o; puts Runner::help(); exit}
end

def parse_known_to(parser, initial_args=ARGV.dup)
    other_args = []                                         # this contains the unknown options
    rec_parse = Proc.new { |arg_list|                       # in_method defined proc
        begin
            parser.parse! arg_list                          # try to parse the arg list
        rescue OptionParser::InvalidOption => e
            other_args += e.args                            # save the unknown arg
            while arg_list[0] && arg_list[0][0] != "-"      # certainly not perfect but
                other_args << arg_list.shift                # quick hack to save any parameters
            end
            rec_parse.call arg_list                         # call itself recursively
        end
    }
    rec_parse.call initial_args                             # start the rec call
    other_args                                              # return the invalid arguments
end

other_options = parse_known_to(opts)
# puts other_options.inspect
# puts options.inspect
# puts ARGV.inspect
command = ARGV.shift
Runner.execute command, options, ARGV

