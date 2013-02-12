#!/usr/bin/env ruby

$:.unshift(File.join(File.dirname(__FILE__), "lib"))

require 'pipette'

puts "running rna_seq pipeline"
pipeline = RnaSeqPipeline.pipeline
# pipeline.default_options
pipeline.parse_input(ARGV)

if !pipeline.options[:output]
  output_prefix = options[:input].split(".")[0..-2].join(".")
  output_prefix = File.join(output_prefix, output_prefix)
  pipeline.options[:output] = output_prefix
end

steps = pipeline.options[:steps] ? pipeline.options[:steps].collect {|step| step.to_sym} : nil

pipeline.run

