#!/usr/bin/env ruby

$:.unshift(File.join(File.dirname(__FILE__), "lib", "vp"))
puts $:.inspect
vcf_file = ARGV[0]
out_file = ARGV[1]
raise "ERROR: usage: variant_filter [VCF_FILE] [OUTPUT_FILE" unless vcf_file && out_file
raise "ERROR: no vcf file found at: #{vcf_file}" unless File.exists? vcf_file

require 'vcf_filter'
require 'benchmark'

time = Benchmark.measure do
  VCFFilter.filter vcf_file, out_file
end
puts time
