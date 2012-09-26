require 'pipette/tools/vcf'

class SnpeffSummarizer
  def run input_filename, options
    raise "ERROR: snpeff txt file not found at: #{input_filename}." unless File.exists? input_filename
    base_name = input_filename.split(".")[0..-2].join(".")
    output_filename = base_name + ".collapse.txt"
    data = read_file(input_filename, options[:raw_vcf_file])

    collapse_file(data, output_filename)

    output_filename
  end

  ALL_DATA = ["Chrom", "Position", "Reference", "Change", "ChangeType", "Homozygous",
              "Quality", "Coverage", "Warnings", "GeneID", "GeneName", "BioType",
              "TranscriptID", "ExonID", "ExonRank", "Effect", "OldAA2NewAA", "OldCondon2NewCondon",
              "CodonNum", "CDS_Size", "CodonsAround", "AAAround", "CustomIntervalID", "Ref_Count", "Change_Count"]

  HEADERS = ["Chrom", "Position", "Reference",
             "Ref_Count", "Change_Count",
             "Homozygous", "Change", "ChangeType",
             "Quality", "Coverage", "GeneID", "GeneName", "BioType",
             "TranscriptID", "Effect", "OldAA2NewAA", "OldCondon2NewCondon"]

  COMBINE = ["TranscriptID", "Effect", "OldAA2NewAA", "OldCondon2NewCondon"]

  def read_file(input_filename, original_vcf_file)
    input_file = File.open(input_filename, 'r')
    puts "Summarizing"
    vcf_hash = {}
    if original_vcf_file
      puts " combining with #{original_vcf_file}"
      vcf = VCF.new(original_vcf_file, {:read_info => false})
      vcf.to_a.each do |vcf_line|
        id = [vcf_line["CHROM"], vcf_line["POS"], vcf_line["REF"], vcf_line["ALT"]].join("_").gsub("chr","")
        vcf_hash[id] = vcf_line
      end
    end
    data = []
    input_file.each_line do |line|
      if line =~ /^#/
        next
      else
        values = line.chomp.split("\t")
        line_data = Hash[ALL_DATA.zip(values)]
        id = [line_data["Chrom"], line_data["Position"], line_data["Reference"], line_data["Change"]].join("_").gsub("chr","")

        if vcf_hash[id]
          ad = vcf_hash[id]["AD"]
          line_data["Ref_Count"] = ad[0]
          line_data["Change_Count"] = ad[1]
        else
          puts "no vcf entry for #{id}"
        end

        data << line_data
      end
    end
    input_file.close
    data
  end

  def collapse_file data, output_filename
    output_file = File.open(output_filename, 'w')
    print_header(HEADERS, output_file)

    previous_id = nil
    previous_data = {}
    data.each do |line_data|
      id = line_data["Chrom"] + ";" + line_data["Position"]
      puts "ERROR: id empty for data" if id.empty?

      if previous_id and (previous_id != id)
        print_data(previous_data, output_file)
        previous_data = {}
      end

      previous_data = copy_data(previous_data, line_data)
      previous_id = id
    end

    print_data(previous_data, output_file)
    output_file.close
  end

  def copy_data old_data, new_data
    data = {}
    HEADERS.each do |header|
      if COMBINE.include? header
        new_data[header] ||= " "
        data[header] = [old_data[header], new_data[header]].compact.join(";")
      else
        data[header] = new_data[header]
      end
    end
    data
  end

  def print_header header, file
    file << header.join("\t") << "\n"
  end

  def print_data data, file
    data = HEADERS.collect {|header| data[header]}
    file << data.join("\t") << "\n"
  end
end
