#!/usr/bin/env ruby

$:.unshift(File.join(File.dirname(__FILE__), "lib"))

require 'optparse'
require 'yaml'
require 'fileutils'
require 'vp'

options = {}
default_options_file = File.expand_path(File.join(File.dirname(__FILE__), "config", "default.yml"))
if File.exists? default_options_file
  options = Hash[YAML::load(open(default_options_file)).map {|k,v| [k.to_sym, v]}]
else
  puts "WARNING: no default configuration found in #{default_options_file}"
end

OptionParser.new do |o|
  o.on('-i', '--input BAM_FILE', 'REQUIRED - Input BAM file to call SNPs on') {|b| options[:input] = b}
  o.on('-r', '--reference FA_FILE', 'REQUIRED - Reference Fasta file for genome') {|b| options[:reference] = b}
  o.on('-o', '--output PREFIX', 'Output prefix to use for generated files') {|b| options[:output] = b}
  o.on('-j', '--cores NUM', "Specify number of cores to run GATK on. Default: #{options[:cores]}") {|b| options[:cores] = b.to_i}
  o.on('-c', '--recalibrate COVARIATE_FILE', "If provided, recalibration will occur using input covariate file. Default: recalibration not performed") {|b| options[:recalibrate] = b}
  o.on('-a', '--annotate GENOME', 'Annotate the SNPs and Indels using Ensembl based on input GENOME. Example Genome: FruitFly') {|b| options[:annotate] = b}
  o.on('--gatk JAR_FILE', String, "Specify GATK installation") {|b| options[:gatk] = b}
  o.on('--snpeff JAR_FILE', String, "Specify snppEff Jar location") {|b| options[:snpeff] = b}
  o.on('--snpeff_config CONFIG_FILE', String, "Specify snppEff config file location") {|b| options[:snpeff_config] = b}
  o.on('--samtools BIN_PATH', String, "Specify location of samtools") {|b| options[:samtools] = b}
  o.on('-q', '--quiet', 'Turn off some output') {|b| options[:verbose] = !b}
  o.on('-s', "--steps STEPS" , Array, 'Specify only which steps of the pipeline should be executed') {|b| options[:steps] = b.collect {|step| step.to_sym} }
  o.on('-y', '--yaml YAML_FILE', String, 'Yaml configuration file that can be used to load options. Command line options will trump yaml options') {|b| options.merge!(Hash[YAML::load(open(b)).map {|k,v| [k.to_sym, v]}]) }
  o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
  o.parse!
end

raise "ERROR - input BAM file required. Use -i parameter, or -h for more info" unless options[:input]
raise "ERROR - reference Fasta file required. Use -r parameter or -h for more info" unless options[:reference]

options[:output] ||= options[:input].split(".")[0..-2].join(".")
options[:samtools] ||= %x[which samtools].chomp

file_options = [:reference, :gatk, :snpeff, :snpeff_config, :samtools]
file_options.each {|opt| options[opt] = File.expand_path(options[opt])}

puts "options used:"
options.each do |option, value|
  puts "#{option} => #{value}"
end

steps = options[:steps] ? options[:steps].collect {|step| step.to_sym} : nil

def check_options options
  raise "ERROR GATK JAR not found at:#{options[:gatk]}." unless File.exists? options[:gatk]
  raise "ERROR samtools not found at:#{options[:samtools]}." unless File.exists? options[:samtools]
  raise "ERROR Reference file not found at:#{options[:reference]}." unless File.exists? options[:reference]
  raise "ERROR Input file not found at:#{options[:input]}" unless File.exists? options[:input]
  if options[:recalibrate]
    raise "ERROR covariate file not found at:#{options[:recalibrate]}." unless File.exists? options[:recalibrate]
  end
end

# ---------------
# MAIN PROGRAM --
# ---------------

check_options options
puts "performing steps: #{steps.join(",")}"

# put analysis in its own sub-directory

puts "creating directory: #{options[:output]}"

FileUtils.mkdir_p prefix unless Dir.exists? options[:output]
options[:output] = File.join(options[:output], options[:output])

pipeline = VariantPipeline.new
pipeline.run options
