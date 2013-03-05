#!/usr/bin/env ruby

require 'yaml'
require 'optparse'
require 'parallel'

options = {}

opts = OptionParser.new do |o|
  o.on('-y', '--yaml YAML_FILE', String, 'Yaml configuration file that can be used to load options. Command line options will trump yaml options') {|b| options.merge!(Hash[YAML::load(open(b)).map {|k,v| [k.to_sym, v]}]) }
  o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
end
opts.parse(ARGV)

puts options.inspect

CUR_DIR = File.expand_path(File.dirname(__FILE__))
PIPELINE = File.expand_path(File.join(CUR_DIR, "..", "rna_seq.rb"))

puts PIPELINE

def assert
  raise "assertion failed !" unless yield
end

class Runner
    HELP = <<-HELP
Usage: rna_seqer [TYPE] [ORDER/STARTING_DIR] <OPTIONS>

TYPE can be one of the following:

  order    - run pipeline on entire order

             expects input to be text file with
             flowcell_id[TAB]flowcell_dir
             for all flowcells in the order

  flowcell - run pipeline on flowcell
             
             expects input to be flowcell_dir
             for flowcell
  

  sample   - run pipeline on a particular sample

             expects input to be full path to 
             sample fastq file

    HELP

  def initialize()
  end

  def load_order(order_file)
    if !File.exists?(order_file)
      puts "ERROR: order file not found"
      puts "#{order_file}"
      exit(1)
    end
    orders = File.open(order_file, 'r').read.split("\n").collect {|l| l.split("\t")}
  end

  def load_sample_report(root_dir)
    sample_report_filename = File.join(root_dir, "Sample_Report.csv")
    if !File.exists?(sample_report_filename)
      puts "ERROR: Sample_Report.csv not found in directory"
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
        samples << sample_hash
      end
    end
    samples
  end

  def join_sample_reports(sample_reports, keys = ["sample name"])
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
      sample['input'] = first_reads
      sample['pair'] = second_reads
      samples << sample
    end
    samples
  end

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

  def load_sample_reports(orders)
    sample_reports = []
    orders.each do |order|
      sample_reports << load_sample_report(order[1])
    end
    samples = join_sample_reports(sample_reports)
    check_samples(samples)
    samples
  end

  def run_order(order_file, args)
    orders = load_order(order_file)
    samples = load_sample_reports(orders)

    run_pipeline(samples, args)
  end

  def run_flowcell(starting_dir, args)
    sample_report = load_sample_report(starting_dir)
    samples = join_sample_reports([sample_report])
    run_pipeline(samples, args)
  end

  def run_pipeline_on_sample(data)
    command = "#{PIPELINE}"
    command += " --name #{data['name']}"
    command += " --input #{data['input']}"
    if data['pair']
      command += " --pair #{data['pair']}"
    end

    command += " #{data['args'].join(" ")}"
    puts command
    puts " ---- "
  end


  def run_pipeline(samples, args)
    processes = 1
    Parallel.each(samples, :in_processes => processes) do |sample|
      input_string = sample['input'].collect{ |i| i['full']}.join(",")
      pair_string = nil
      if sample['pair'] and !sample['pair'].empty?
        pair_string = sample['pair'].collect{ |i| i['full']}.join(",")
      end

      sample_data = {'name' => sample['name'], 'input' => input_string, 'pair' => pair_string, 'args' => args}
      run_pipeline_on_sample(sample_data)
    end
  end

  def self.execute(command, args)
    case command
    when "order"
      self.new().run_order(args.shift, args)
    when "flowcell"
      self.new().run_flowcell(args.shift, args)
    when "sample"
    else
      puts Runner::HELP
    end
  end
end

command = ARGV.shift
Runner.execute command, ARGV

