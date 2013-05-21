#!/usr/bin/env ruby

$LOAD_PATH.unshift(File.dirname(__FILE__))

require 'sample_report'

bam_dir = ARGV[0]

output_dir = File.join(bam_dir, "joined_bams")
system("mkdir -p #{output_dir}")

sample_report_filename = File.join(bam_dir, "Sample_Report.csv")
if !File.exists? sample_report_filename
  puts "ERROR: no sample report in #{path}"
  exit(1)
end

sample_report = SampleReport.new(sample_report_filename)

bam_filenames = Dir.glob(File.join(bam_dir, "**", "*.group.reorder.bam"))
# bam_filenames = Dir.glob(File.join(bam_dir, "**", "*.rmdup.group.bam"))

puts "#{bam_filenames.size} found"

samples = {}

bam_filenames.each do |bam_filename|
  fastq_filename = File.basename(File.dirname(bam_filename)) + ".fastq.gz"
  puts fastq_filename

  fastq_data = sample_report.data_for(fastq_filename)
  sample_name = fastq_data["sample name"]
  puts sample_name

  samples[sample_name] ||= []
  samples[sample_name] << bam_filename

end

puts "#{samples.keys.size} samples found"

samples.each do |name, files|
  output_filename = File.join(output_dir, name + ".bam")

  command = "samtools merge #{output_filename} #{files.join(" ")}"
  puts command
  system command

  command = "samtools index #{output_filename}"
  system command
end


