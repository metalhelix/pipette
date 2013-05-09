#!/usr/bin/env ruby

CUR_DIR = File.dirname(__FILE__)
$LOAD_PATH.unshift(CUR_DIR)

require 'sample_report'

require 'parallel'

PIPETTE = File.join(CUR_DIR, "..", "pe_bwa.rb")

if !File.exists? PIPETTE
  puts "ERROR: cannot find pipette at:"
  puts " #{PIPETTE}"
  exit(1)
end

WORKING_DIR = Dir.pwd
PROCESSES = 2


fastq_dir = File.expand_path(ARGV[0])
config_file = ARGV[1]

if !fastq_dir or !File.exists?(fastq_dir) or !config_file or !File.exists?(config_file)
  puts "ERROR: invalid input"
  puts "Please run with the fastq directory and the pipette config file."
  puts "Example: "
  puts "  ./pipette/bulk/bwa_align.rb /n/analysis/Mak/huy/MOLNG-395/C1TN6ACXX/ pe_bwa_config.yml"
  exit(1)
end



output_dir = File.join(WORKING_DIR, "align")
system("mkdir -p #{output_dir}")


# reject second read as they will be handled below
sequence_filenames = Dir.glob(File.join(fastq_dir, "s_*_1_*.fastq.gz"))
sequence_filenames.reject! {|name| name =~ /Undetermined/}

puts "Found #{sequence_filenames.size} sequences"


sample_report_filename = File.join(fastq_dir, "Sample_Report.csv")
if !File.exists? sample_report_filename
  puts "ERROR: no sample report in #{path}"
  exit(1)
end
sample_report = SampleReport.new(sample_report_filename)

Parallel.each(sequence_filenames, :in_processes => PROCESSES) do |sequence_filename|
  sequence_basename = File.basename(sequence_filename)
  # original_name = sequence_basename.gsub(".trim","")
  # sequence_data = report.data_for original_name

  # name = sequence_data["sample name"]
  # name = name.downcase.gsub(" ", "")
  filename = File.basename(sequence_basename, File.extname(sequence_basename)).split(".")[0]
  input_filename = sequence_filename
  #TODO fix badness
  pair_filename = filename.split("_")[0..1].join("_") + "_2_" + filename.split("_")[-1] + ".fastq.gz"
  pair_filename = File.join(File.dirname(sequence_filename), pair_filename)
  # pair = sequence_filename.split("_")[0] + "_2.fastq.gz"
  if !File.exists?(pair_filename)
    puts "ERROR: file not found #{pair_filename}"
    next
  end


  # fastq_data = sample_report.data_for(pair_filename)
  # sample_name = fastq_data["sample name"]
  # puts sample_name
  sample_name = filename

  full_output_dir = File.join(output_dir, "#{sample_name}", "#{sample_name}")

  if File.exists?(File.dirname(full_output_dir))
    puts "skipping #{full_output_dir}"
    next
  end

  # puts sample_name
  # puts input_filename
  # puts pair_filename

  command = "#{PIPETTE} -y #{config_file} --input #{input_filename} --pair #{pair_filename} --output #{full_output_dir} --name #{sample_name}"
  puts command
  result = %x[#{command}]
  puts result
end


