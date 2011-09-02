#!/usr/bin/env ruby

$:.unshift(File.join(File.dirname(__FILE__), "lib"))

require 'vp'

puts "running variant pipeline"
pipeline = VariantPipeline.pipeline
pipeline.default_options
pipeline.parse_input(ARGV)

if !pipeline.options[:output]
  output_prefix = options[:input].split(".")[0..-2].join(".")
  output_prefix = File.join(output_prefix, output_prefix)
  pipeline.options[:output] = output_prefix
end

pipeline.options[:samtools] ||= %x[which samtools].chomp

file_options = [:reference, :gatk, :snpeff, :snpeff_config, :samtools]
file_options.each {|opt| pipeline.options[opt] = File.expand_path(options[opt])}

puts "options used:"
pipeline.options.each do |option, value|
  puts "#{option} => #{value}"
end

steps = pipeline.options[:steps] ? pipeline.options[:steps].collect {|step| step.to_sym} : nil

# ---------------
# MAIN PROGRAM --
# ---------------

puts "performing steps: #{steps.join(",")}"

# put analysis in its own sub-directory
output_dir = File.dirname(pipeline.options[:output])
puts "creating directory: #{output_dir}"

FileUtils.mkdir_p output_dir unless Dir.exists? output_dir

pipeline.run

