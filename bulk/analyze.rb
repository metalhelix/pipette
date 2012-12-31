#!/usr/bin/env ruby

# recommended running:
# nohup ./analyze.rb /path/to/align 2>&1 > analyze.out.log &


input_bam_dir = ARGV[0]

if !input_bam_dir
  puts "ERROR: run using ./analyze.rb [path/to/align/dir]"
  exit(1)
end

input_bam_dir = File.expand_path(input_bam_dir)

cur_dir = File.expand_path(File.dirname(__FILE__))
# filter = "*.filter.bam"
filter = "*.bam"

bam_files = Dir.glob(File.join(input_bam_dir, "**", filter))

puts "Found #{bam_files.size} matching bam files"

PIPETTE = File.join(cur_dir, 'pipette', 'vp.rb')
CONFIG_FILE = File.join(cur_dir, "variant_config.yml")

if !File.exists?(PIPETTE)
  puts "ERROR: can't find file #{PIPETTE}"
  exit(1)
end

if !File.exists?(CONFIG_FILE)
  puts "ERROR: can't find file #{CONFIG_FILE}"
  exit(1)
end

bam_files.each do |input_bam_file|

  prefix = File.basename(input_bam_file).split(".")[0]
  puts prefix

  output = File.join(cur_dir, "analyze", prefix, prefix)

  if File.exists?(File.dirname(output))
    puts "EXISTS: #{output}"
    next
  end

  puts "Creating output here:#{output}"

  command = "#{PIPETTE} -y #{CONFIG_FILE} --input #{input_bam_file} --output #{output}"
  puts command
  result = %x[#{command}]
  puts result
end

