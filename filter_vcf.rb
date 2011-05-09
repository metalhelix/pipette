#!/usr/bin/env ruby

$:.unshift(File.dirname(__FILE__))

vcf_file = ARGV[0]
out_file = ARGV[1]
raise "ERROR: usage: variant_filter [VCF_FILE] [OUTPUT_FILE" unless vcf_file && out_file
raise "ERROR: no vcf file found at: #{vcf_file}" unless File.exists? vcf_file

require 'vcf'
require 'benchmark'

time = Benchmark.measure do
vcf = VCF.new vcf_file

vcf.copy_if(out_file) do |vcf_hash| 
  rtn = vcf_hash["DP2"] > 10 
  if rtn && vcf_hash["genotype"] == "isHet"
    val  = (vcf_hash["AD"][1].to_f / vcf_hash["AD"][0].to_f)
    rtn = val > 2.0
    #puts "low divide: #{val}" unless rtn
  end
  rtn
end

vcf.close
end
puts time
