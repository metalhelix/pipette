#!/usr/bin/env ruby

$:.unshift(File.join(File.dirname(__FILE__), "lib"))


require 'optparse'
require 'yaml'
require 'fileutils'
require 'vp'

@valid_steps = [:realign, :recalibrate, :call, :filter, :annotate]
$options = {}
$options[:recalibrate] = false
$options[:verbose] = true
$options[:cores] = 4
$options[:gatk] = "/n/site/inst/Linux-x86_64/bioinfo/GATK/GenomeAnalysisTK-1.0.5315/GenomeAnalysisTK.jar"
$options[:steps] = @valid_steps

OptionParser.new do |o|
  o.on('-i', '--input BAM_FILE', 'REQUIRED - Input BAM file to call SNPs on') {|b| $options[:input] = b}
  o.on('-r', '--reference FA_FILE', 'REQUIRED - Reference Fasta file for genome') {|b| $options[:reference] = b}
  o.on('-o', '--output PREFIX', 'Output prefix to use for generated files') {|b| $options[:output] = b}
  o.on('-j', '--cores NUM', "Specify number of cores to run GATK on. Default: #{$options[:cores]}") {|b| $options[:cores] = b.to_i}
  o.on('-c', '--recalibrate COVARIATE_FILE', "If provided, recalibration will occur using input covariate file. Default: recalibration not performed") {|b| $options[:recalibrate] = b}
  o.on('-a', '--annotate GENOME', 'Annotate the SNPs and Indels using Ensembl based on input GENOME. Example Genome: FruitFly') {|b| $options[:annotate] = b}
  o.on('-g', '--gatk JAR_FILE', "Specify GATK installation") {|b| $options[:gatk] = b}
  o.on('-q', '--quiet', 'Turn off some output') {|b| $options[:verbose] = !b}
  o.on('-s', "--steps #{@valid_steps.join(",")}", Array 'Specify only which steps of the pipeline should be executed') {|b| $options[:steps] = b.collect {|step| step.to_sym} }
  o.on('-y', '--yaml YAML_FILE', String, 'Yaml configuration file that can be used to load options. Command line options will trump yaml options') {|b| $options.merge!(Hash[YAML::load(open(b)).map {|k,v| [k.to_sym, v]}]) }
  o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
  o.parse!
end

raise "ERROR - input BAM file required. Use -i parameter, or -h for more info" unless $options[:input]
raise "ERROR - reference Fasta file required. Use -r parameter or -h for more info" unless $options[:reference]

$options[:log_dir] ||= "log"
$options[:out_dir] ||= "out"
$options[:output] ||= $options[:input].split(".")[0..-2].join(".")


puts "options used:\n#{$options.inspect}"

@steps = $options[:steps] ? $options[:steps].collect {|step| step.to_sym} : nil

def check_options
  raise "ERROR GATK JAR not found at:#{$options[:gatk]}." unless File.exists? $options[:gatk]
  raise "ERROR Reference file not found at:#{$options[:reference]}." unless File.exists? $options[:reference]
  raise "ERROR Input file not found at:#{$options[:input]}" unless File.exists? $options[:input]
  if $options[:recalibrate]
    raise "ERROR covariate file not found at:#{$options[:recalibrate]}." unless File.exists? $options[:recalibrate]
  end

end

def check_steps
  if !@steps || @steps.empty?
    raise "ERROR no steps provided"
  end
  @steps.each do |step|
    unless @valid_steps.include? step
      raise "ERROR: #{step} not a valid step.\nvalid steps: #{@valid_steps.join(",")}"
    end
  end
end

# ---------------
# MAIN PROGRAM --
# ---------------

check_options
check_steps
puts "performing steps: #{@steps.join(",")}"

# put analysis in its own sub-directory
prefix = $options[:output]
FileUtils.mkdir_p prefix unless Dir.exists? prefix
prefix = "#{prefix}/#{prefix}"

pipeline = Pipeline.new $options
pipeline.prefix = prefix
pipeline.steps = @steps

variant_bam_file = pipeline.realign $options[:input], prefix

should_recal = $options[:recalibrate]
if should_recal
  variant_bam_file = pipeline.recalibrate $options[:input], variant_bam_file, prefix
end

vcf_files = []
snps_vcf_file = pipeline.call_snps variant_bam_file, prefix

vcf_files << snps_vcf_file

indels_vcf_file = pipeline.call_indels variant_bam_file, prefix

vcf_files << indels_vcf_file

filtered_vcf_files = pipeline.filter vcf_files

if $options[:annotate]
  puts "annoting: #{filtered_vcf_files.inspect}"
  pipeline.annotate filtered_vcf_files
end

puts "Pipeline complete!"

