#!/usr/bin/env ruby

$:.unshift(File.join(File.dirname(__FILE__), "lib"))

require 'pipette'

ALL_PIPELINES = { "variant" => {:name => "variant pipeline", :class => VariantPipeline},
                  "pe_bwa" => {:name => "paired-end bwa pipeline", :class => BwaPipeline}
                }

input_pipeline_name = ARGV.shift

if !input_pipeline_name
  puts "ERROR: No pipeline name provided"
  puts "please run using ./pipette.rb <pipeline_name>"
  puts "Valid pipeline names are:"
  ALL_PIPELINES.keys.each {|k| puts "\t#{k}"}
  exit
end

pipeline_data = ALL_PIPELINES[input_pipeline_name]

if !pipeline_data
  puts "#{input_pipeline_name} not a valid pipeline"
  puts "Valid pipeline names are:"
  ALL_PIPELINES.keys.each {|k| puts "\t#{k}"}
  exit
end

puts "running #{pipeline_data[:name]}"

pipeline = pipeline_data[:class].pipeline

pipeline.default_options
pipeline.parse_input(ARGV)

if !pipeline.options[:output]
  output_prefix = options[:input].split(".")[0..-2].join(".")
  output_prefix = File.join(output_prefix, output_prefix)
  pipeline.options[:output] = output_prefix
end

steps = pipeline.options[:steps] ? pipeline.options[:steps].collect {|step| step.to_sym} : nil

pipeline.run


