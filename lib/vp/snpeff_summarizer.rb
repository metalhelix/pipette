
class SnpeffSummarizer
  def initialize(options = {})
  end

  def run input_filename, options
    raise "ERROR: snpeff txt file not found at: #{input_filename}." unless File.exists? input_filename
    base_name = input_filename.split(".")[0..-2].join(".")
    output_filename = base_name + ".collapse.txt"

    collapse_file(input_filename, output_filename)

    output_filename
  end

  ALL_DATA = ["Chrom", "Position", "Reference", "Change", "ChangeType", "Homozygous",
              "Quality", "Coverage", "Warnings", "GeneID", "GeneName", "BioType",
              "TranscriptID", "ExonID", "ExonRank", "Effect", "OldAA2NewAA", "OldCondon2NewCondon",
              "CodonNum", "CDS_Size", "CodonsAround", "AAAround", "CustomIntervalID"]

  HEADERS = ["Chrom", "Position", "Reference", "Change", "ChangeType",
             "Quality", "Coverage", "GeneID", "GeneName", "BioType",
             "TranscriptID", "Effect", "OldAA2NewAA", "OldCondon2NewCondon"]

  COMBINE = ["TranscriptID", "Effect", "OldAA2NewAA", "OldCondon2NewCondon"]

  def collapse_file input_filename, output_filename
    output_file = File.open(output_filename, 'w')
    print_header(HEADERS, output_file)
    input_file = File.open(input_filename, 'r')
    previous_id = nil
    previous_data = {}
    input_file.each_line do |line|
      if line =~ /^#/
        previous_id = nil
      else
        values = line.chomp.split("\t")
        line_data = Hash[ALL_DATA.zip(values)]
        id = line_data["Chrom"] + ";" + line_data["Position"]
        puts "ERROR: id empty for data" if id.empty?

        if previous_id and (previous_id != id)
          print_data(previous_data, output_file)
          previous_data = {}
        end

        previous_data = copy_data(previous_data, line_data)
        previous_id = id
      end
    end

    print_data(previous_data, output_file)
    output_file.close
    input_file.close
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
