#!/usr/bin/env ruby

$:.unshift(File.dirname(__FILE__))

require 'gatk'
require 'optparse'

$options = {}
$options[:recal] = false
$options[:verbose] = true
$options[:cores] = 4

OptionParser.new do |o|
  o.on('-i', '--input VCF_FILE', 'REQUIRED - Input VCF file to filter') {|b| $options[:input] = b}
  o.on('-r', '--reference FA_FILE', 'REQUIRED - Reference Fasta file for genome') {|b| $options[:reference] = b}
  o.on('-o', '--output PREFIX', 'Output prefix to use for generated files') {|b| $options[:output] = b}
  o.on('-g', '--gatk JAR_FILE', "Specify GATK installation. Default: #{$options[:jar]}") {|b| $options[:jar] = b}
  o.on('-q', '--quiet', 'Turn off some output') {|b| $options[:verbose] = !b}
  o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
  o.parse!
end

raise "ERROR - input vcf file required. Use -i parameter or -h for more info" unless $options[:input]
raise "ERROR - reference Fasta file required. Use -r parameter or -h for more info" unless $options[:reference]


$options[:log_dir] ||= "log"
$options[:out_dir] ||= "out"
$options[:output] ||= $options[:input].split(".")[0]
prefix = $options[:output]

def report status
  puts "#{Time.now} - " + status if $options[:verbose]
end

def valid_genome? genome
  valid_genomes = ["FruitFly"]

  valid_genomes.include? genome
end

def check_options
  raise "ERROR Reference file not found at #{$options[:reference]}" unless File.exists? $options[:reference]
  raise "ERROR Input file not found at #{$options[:input]}" unless File.exists? $options[:input]
  if $options[:recal]
    raise "ERROR covariate file not found at #{$options[:recal]}" unless File.exists? $options[:recal]
  end
  if $options[:annotate]
    raise "ERROR invalid genome provided: #{$options[:annotate]}" unless valid_genome? $options[:annotate]

  end
end

def database_for genome
  case genome
  when /[Ff]ly/
    "drosophila_melanogaster_core_57_513b"
  else
    raise "ERROR invalid genome provided"
  end
end

check_options

gatk = GATK.new $options
report "Starting filtering"
out_file = "#{prefix}.filtered.vcf"
params = {"-T" => "VariantFiltration", 
          "-B:variant,VCF" => $options[:input], 
          "-o" => out_file,
          "--filterName" => "FilterDP",
          "--filterExpression" => "\"DP <= 10\""}
gatk.execute params
