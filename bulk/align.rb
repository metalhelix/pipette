#!/usr/bin/env ruby

$LOAD_PATH.unshift(File.dirname(__FILE__))

require 'sample_report'

require 'parallel'

CUR_DIR = File.expand_path(File.dirname(__FILE__))
$:.unshift(CUR_DIR)

PROCESSES = 2


fastq_dir = File.expand_path(ARGV[0])
output_dir = File.join(CUR_DIR, "align")
system("mkdir -p #{output_dir}")

PIPETTE = File.join(CUR_DIR, "..", "pe_bwa.rb")
CONFIG_FILE = File.join(CUR_DIR, "pe_bwa_config.yml")


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

  command = "#{PIPETTE} -y #{CONFIG_FILE} --input #{input_filename} --pair #{pair_filename} --output #{full_output_dir} --name #{sample_name}"
  puts command
  result = %x[#{command}]
  puts result
end


