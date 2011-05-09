#!/usr/bin/env ruby

$:.unshift(File.dirname(__FILE__))

require 'gatk'
require 'optparse'

$options = {}
$options[:recal] = false
$options[:verbose] = true
$options[:cores] = 4

OptionParser.new do |o|
  o.on('-i', '--input BAM_FILE', 'REQUIRED - Input BAM file to call SNPs on') {|b| $options[:input] = b}
  o.on('-r', '--reference FA_FILE', 'REQUIRED - Reference Fasta file for genome') {|b| $options[:reference] = b}
  o.on('-o', '--output PREFIX', 'Output prefix to use for generated files') {|b| $options[:output] = b}
  o.on('-j', '--cores NUM', "Specify number of cores to run GATK on. Default: #{$options[:cores]}") {|b| $options[:cores] = b.to_i}
  o.on('-c', '--recalibrate COVARIATE_FILE', "If provided, recalibration will occur using input covariate file. Default: recalibration not performed") {|b| $options[:recal] = b}
  o.on('-a', '--annotate GENOME', 'Annotate the SNPs and Indels using Ensembl based on input GENOME. Example Genome: FruitFly') {|b| $options[:annotate] = b}
  o.on('-g', '--gatk JAR_FILE', "Specify GATK installation. Default: #{$options[:jar]}") {|b| $options[:jar] = b}
  o.on('-q', '--quiet', 'Turn off some output') {|b| $options[:verbose] = !b}
  o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
  o.parse!
end

raise "ERROR - input BAM file required. Use -i parameter, or -h for more info" unless $options[:input]
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
  raise "ERROR GATK not found at #{$options[:gatk]}" unless File.exists? $options[:gatk]
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

def execute command, fork = false
  report command
  #result = %x[#{command}]
  result = system(command)
  #report result
end


check_options

gatk = GATK.new $options

annotate_script = nil
if $options[:annotate]
  annotate_script = File.join(File.dirname(__FILE__), "annotate", "annotateStrainSNVDiffs.pl")
  raise "ERROR: annotation script not found at: #{annotate_script}" unless File.exists? annotate_script
end


report "Starting realignment interval creation"
intervals_file = "#{prefix}.intervals"
params = {"-T" => "RealignerTargetCreator", 
          "-I" => $options[:input], 
          "-o" => intervals_file}

gatk.execute params


report "Starting realignment using intervals"
realign_bam_file = "#{prefix}.realigned.bam"
params = {"-T" => "IndelRealigner", 
          "-I" => $options[:input], 
          "-targetIntervals" => intervals_file, 
          "-o" => realign_bam_file}

gatk.execute params

report "Starting index of realigned BAM with Samtools"
command = "samtools index #{realign_bam_file}"
execute command

variant_bam_file = realign_bam_file

should_recal = $options[:recal]
if should_recal 
  report "Starting recalibration"
  covar_table_file = "#{prefix}.covariate_table.csv"

  params = {"-T" => "CountCovariates",
            "-I" => $options[:input],
            "-recalFile" => covar_table_file,
            "-cov" => ["ReadGroupcovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate"]}
  gatk.execute params

  recal_bam_file = "#{prefix}.realigned.recal.bam"

  params = {"-T" => "TableRecalibration",
            "-I" => realign_bam_file,
            "-recalFile" => covar_table_file,
            "--out" => recal_bam_file}
  gatk.execute params

  cleaned_bam_file = "#{prefix}.realigned.recal.sorted.bam"
  report "Starting sort and index recalibrated BAM"
  command = "samtools sort #{recal_bam_file} #{cleaned_bam_file}"
  execute command
  command = "samtools index #{cleaned_bam_file}"
  execute command

  variant_bam_file = cleaned_bam_file
end

report "Starting SNP calling"
snps_vcf_file = "#{prefix}.snps.vcf"

params = {"-T" => "UnifiedGenotyper",
          "-I" => variant_bam_file,
          "-o" => snps_vcf_file,
          "-stand_call_conf" => "30.0",
          "-stand_emit_conf" => "10.0"}
gatk.execute params

report "Starting Indel calling"
indels_vcf_file = "#{prefix}.indels.vcf"

params = {"-T" => "UnifiedGenotyper",
          "-I" => variant_bam_file,
          "-o" => indels_vcf_file,
          "-glm" => "DINDEL",
          "-stand_call_conf" => "30.0",
          "-stand_emit_conf" => "10.0"}
gatk.execute params

should_annotate = $options[:annotate]
if should_annotate

  genome = $options[:annotate]
  database = database_for genome
  port = "5306"
  site = "ensembldb.ensembl.org"

  report "Starting SNP annotation"
  snp_annotation_file = "#{prefix}.snp.annotate.csv"
  snp_annotation_log_file = File.join(File.dirname(__FILE__),"log","#{prefix}.snp.annotate.log")

  base_command = "#{annotate_script} "
  database_options = " #{genome} #{database} #{port} #{site}"
  vcf_file = snps_vcf_file 
  command = base_command + vcf_file + database_options + " 1> #{snp_annotation_file} 2> #{snp_annotation_log_file}"
  execute command, true

  report "Starting Indel annotation"
  indel_annotation_file = "#{prefix}.indel.annotate.csv"
  indel_annotation_log_file = "#{prefix}.indel.annotate.log"
  vcf_file = indels_vcf_file
  command = base_command + vcf_file + database_options + " 1> #{indel_annotation_file} 2> #{indel_annotation_log_file}" 
  execute command, true
end

Process.waitall
report "Pipeline complete!"

