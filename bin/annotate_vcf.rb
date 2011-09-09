#!/usr/bin/env ruby

$:.unshift(File.join(File.dirname(__FILE__), "..", "lib"))

require 'pipette/tools/annotator'
require 'benchmark'
require 'optparse'
require 'yaml'

$options = {}
$options[:snpeff] = "/n/site/inst/Linux-x86_64/bioinfo/snpEff/current/snpEff.jar"
$options[:snpeff_config] = "/n/site/inst/Linux-x86_64/bioinfo/snpEff/current/snpEff.config"
$options[:annotate] = "dm5.34"

OptionParser.new do |o|
  o.on('-i', '--input VCF_FILE', 'REQUIRED - Input VCF file to annotate SNPs') {|b| $options[:input] = b}
  o.on('-y', '--yaml YAML_FILE', String, 'Yaml configuration file that can be used to load options. Command line options will trump yaml options') {|b| $options.merge!(Hash[YAML::load(open(b)).map {|k,v| [k.to_sym, v]}]) }
  o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
  o.parse!
end


time = Benchmark.measure do
  annotator = Annotator.new
  output_file = annotator.run($options[:input], $options)
end
puts time
