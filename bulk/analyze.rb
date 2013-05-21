#!/usr/bin/env ruby

# recommended running:
# nohup ./analyze.rb /path/to/align 2>&1 > analyze.out.log &


input_bam_dir = ARGV[0]
config_file = ARGV[1]

if !input_bam_dir
  puts "ERROR: run using ./analyze.rb [path/to/align/dir]"
  exit(1)
end

input_bam_dir = File.expand_path(input_bam_dir)

cur_dir = File.expand_path(File.dirname(__FILE__))
# filter = "*.filter.bam"
filter = "*.bam"

output_dir = File.expand_path(Dir.getwd)

bam_files = Dir.glob(File.join(input_bam_dir, "**", filter))

puts "Found #{bam_files.size} matching bam files"

PIPETTE = File.join(cur_dir, "..",  'vp.rb')

if !File.exists?(PIPETTE)
  puts "ERROR: can't find file #{PIPETTE}"
  exit(1)
end

if !File.exists?(config_file)
  puts "ERROR: can't find file #{config_file}"
  exit(1)
end

bam_files.each do |input_bam_file|

  prefix = File.basename(input_bam_file).split(".")[0]
  puts prefix

  output = File.join(output_dir, "analyze", prefix, prefix)

  if File.exists?(File.dirname(output))
    puts "EXISTS: #{output}"
    next
  end

  puts "Creating output here:#{output}"

  command = "#{PIPETTE} -y #{config_file} --input #{input_bam_file} --output #{output}"
  puts command
  result = %x[#{command}]
  puts result
end

