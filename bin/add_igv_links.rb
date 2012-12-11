#!/usr/bin/env ruby

input_snpeff_txt_file = ARGV[0]

output_filename = File.join(File.dirname(input_snpeff_txt_file), File.basename(input_snpeff_txt_file, File.extname(input_snpeff_txt_file)) + ".links.xls")

output_file = File.open(output_filename, 'w')

File.open(input_snpeff_txt_file, 'r') do |file|
  header = nil
  file.each_line do |line|
    fields = line.chomp.split("\t")
    if !header
      fields << "link"
      header = fields
    else
      while fields.size < (header.size - 1)
        fields << " "
      end
      if fields[6][0..0] == "+" or fields[6][0..0] == "-"
        fields[6] = "_ #{fields[6]}"
      end
      loc = "#{fields[0]}:#{fields[1]}"
      link = "=HYPERLINK(\"http://localhost:60151/goto?locus=#{loc}\", \"#{loc}\")"
      fields << link
    end
    output_file.puts fields.join("\t")
  end
end


output_file.close
