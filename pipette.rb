#!/usr/bin/env ruby

$:.unshift(File.join(File.dirname(__FILE__), "lib"))

require 'pipette'

# Get the first argument. Shift so
# pipelines don't receive it
input_pipeline_name = ARGV.shift

#ALL_PIPELINES is defined in lib/pipette/pipelines

# Fail if no argument is provided
if !input_pipeline_name
  puts "ERROR: No pipeline name provided"
  puts "please run using ./pipette.rb <pipeline_name>"
  puts "Valid pipeline names are:"
  ALL_PIPELINES.keys.each {|k| puts "\t#{k}"}
  exit
end

# Acquire pipeline data for that name
pipeline_data = ALL_PIPELINES[input_pipeline_name]

# Fail if not a valid pipeline name
if !pipeline_data
  puts "#{input_pipeline_name} not a valid pipeline"
  puts "Valid pipeline names are:"
  ALL_PIPELINES.keys.each {|k| puts "\t#{k}"}
  exit
end

puts "running #{pipeline_data[:name]}"

# Pull out the pipeline from the pipeline class
pipeline = pipeline_data[:class].pipeline

# Rest of the code below should probably be moved
# to inside pipette framework
pipeline.default_options
pipeline.parse_input(ARGV)

if !pipeline.options[:output]
  output_prefix = options[:input].split(".")[0..-2].join(".")
  output_prefix = File.join(output_prefix, output_prefix)
  pipeline.options[:output] = output_prefix
end

steps = pipeline.options[:steps] ? pipeline.options[:steps].collect {|step| step.to_sym} : nil

pipeline.run


