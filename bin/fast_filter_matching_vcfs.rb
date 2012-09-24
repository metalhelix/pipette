#!/usr/bin/env ruby

require 'set'
require 'optparse'

snps_1_filename = ARGV[0]
snps_2_filename = ARGV[1]

if !snps_1_filename or !snps_2_filename
  puts "ERROR: need two vcf files"
  raise "invalid inputs"
end

options = {:shared => false, :out => File.basename(snps_1_filename)}
opts = OptionParser.new do |o|
  o.banner = "Usage: fast_filter_matching_vcfs path/to/mutant/snps path/to/control/snps [options]"
  o.on('-o', '--out /file/path', String, 'specify the output path') {|b| options[:out] = File.expand_path(b)}
  o.on('-s', '--shared', 'return snps shared instead of unique snps') {|b| options[:shared] = true}
  # o.on('-y', '--yaml YAML_FILE', String, "Yaml configuration file that can be used to load options.","Command line options will trump yaml options") {|b| options.merge!(Hash[YAML::load(open(b)).map {|k,v| [k.to_sym, v]}]) }
  o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
end
opts.parse!
@options = options

puts @options.inspect

shared = @options[:shared]

snp_2_name = File.basename(snps_2_filename).split(".")[0..-2].join(".")

qualifier = shared ? ".shared_with_" : ".unique_from_"

output_name = File.join(@options[:out], File.basename(snps_1_filename).split(".")[0..-2].join(".") + "#{qualifier}#{snp_2_name}.vcf")
puts output_name

def parse_vcf_file filename
  contents = File.open(filename, 'r').read.split("\n")
  while(contents[0][0..1] == "##")
    contents.shift
  end
  header = contents.shift.split("\t")
  data = contents.collect {|line| Hash[header.zip(line.split("\t"))]}
  data
end

def hash_snp snp
  hash_components = []
  ["#CHROM", "POS", "ALT"].each do |field|
    hash_components << snp[field]
  end
  hash_components.join("_")
end

def add_header output_file, input_filename
  contents = File.open(input_filename, 'r').read.split("\n")
  while(contents[0][0..1] == "##")
    head = contents.shift
    output_file.puts head
  end
end

def make_set snps_hash
  snp_hashes = snps_hash.collect {|sh| hash_snp(sh)}
  snp_set = Set.new(snp_hashes)
end

puts "parsing #{snps_1_filename}"
snps_1 = parse_vcf_file snps_1_filename
puts "parsing #{snps_2_filename}"
snps_2 = parse_vcf_file snps_2_filename

puts "hashing #{snps_2_filename}"
snps_2_set = make_set snps_2

puts "excluding variants in #{snps_1_filename} that are found in #{snps_2_filename}"

output = File.open(output_name, 'w')

add_header output, snps_1_filename
output.puts(snps_1.first.keys.join("\t"))

snps_1.each do |snp|
  included = (snps_2_set.include?(hash_snp(snp)))
  keep = shared ? included : !included
  if keep
    output.puts(snp.values.join("\t"))
  end
end

output.close


