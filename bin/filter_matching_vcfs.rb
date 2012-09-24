#!/usr/bin/env ruby

snps_1_filename = ARGV[0]
snps_2_filename = ARGV[1]


if !snps_1_filename or !snps_2_filename
  puts "ERROR: need two vcf files"
  raise "invalid inputs"
end

output_name = snps_1_filename.split(".")[0..-2].join(".") + ".unique.vcf"

def parse_vcf_file filename
  contents = File.open(filename, 'r').read.split("\n")
  while(contents[0][0..1] == "##")
    contents.shift
  end
  header = contents.shift.split("\t")
  data = contents.collect {|line| Hash[header.zip(line.split("\t"))]}
  data
end

def same_snp snp_1, snp_2
  same = true
  ["#CHROM", "POS", "ALT"].each do |field|
    if !snp_1[field]
      puts "ERROR: missing #{field}"
      puts snp_1.inspect
    end
    if !snp_2[field]
      puts "ERROR: missing #{field}"
      puts snp_2.inspect
    end
    if snp_1[field] != snp_2[field]
      same = false
      return same
    end
  end
  same
end

def add_header output_file, input_filename
  contents = File.open(input_filename, 'r').read.split("\n")
  while(contents[0][0..1] == "##")
    head = contents.shift
    output_file.puts head
  end
end

snps_1 = parse_vcf_file snps_1_filename
snps_2 = parse_vcf_file snps_2_filename

puts "excluding variants in #{snps_1_filename} that are found in #{snps_2_filename}"

output = File.open(output_name, 'w')

add_header output, snps_1_filename
output.puts(snps_1.first.keys.join("\t"))

snps_1.each do |snp|
  keep = true
  snps_2.each do |snp2|
    if same_snp(snp,snp2)
      keep = false
      break
    end
  end

  if keep
    output.puts(snp.values.join("\t"))
  end
end

output.close


