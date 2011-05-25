#!/usr/bin/env ruby

require "./variant_pipeline/lib/vp/vcf.rb"

csv_filename = ARGV[0]
vcf_filename = ARGV[1]
output_filename = "out.vcf"
outputfile = File.new(output_filename, 'w')

raise "ERROR: csv file missing" unless csv_filename
raise "ERROR: vcf file missing" unless vcf_filename

csv_lines = File.open(csv_filename, 'r') {|csv_file| csv_lines = csv_file.readlines}

header = csv_lines.shift.chomp.gsub(/"/,"").split(",")
csv_data = csv_lines.collect do |line|
  data = Hash[header.zip(line.chomp.gsub(/"/,"").split(","))]
  data
end

csv_data_copy = csv_data

vcf_file = VCF.new vcf_filename

vcf_file.each do |vcf|
  csv_data.each_with_index do |csv,index|
    csv_id = "#{csv["Chrom"]};#{csv["Position"]}"
    vcf_id = "#{vcf["CHROM"]};#{vcf["POS"]}"
    if csv_id == vcf_id
      puts "Match found #{csv_id} - #{vcf_id}"
      outputfile << vcf.inspect << "\n"
      csv_data_copy.delete csv 
    end
    if csv_data_copy.empty?
      break
    end
  end
end

puts "# of csv left: #{csv_data_copy.size}"
outputfile.close
