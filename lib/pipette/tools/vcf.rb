
class VCF
  def initialize filename
    filename = File.expand_path(filename)
    raise "VCF File not found at: #{filename}" unless File.exists? filename
    @file = File.open(filename,'r')
  end

  def to_h header, line
    geno_data_key = header[header.index("FORMAT") + 1]
    values = header.zip(line.chomp.split("\t"))
    hash_values = Hash[*values.flatten]
    formats = hash_values["FORMAT"].split(":")
    values = hash_values[geno_data_key].split(":")
    formats = formats.zip values
    hash_formats = Hash[*formats.flatten]
    hash_values.merge! hash_formats

    hash_values["genotype"] = case hash_values["GT"]
    when "0/0" then "isHomRef"
    when "0/1" then "isHet"
    when "1/1" then "isHomVar"
    else "unknown"
    end

    hash_values["AD"] = hash_values["AD"].split(",").collect {|a| a.to_i}
    hash_values["PL"] = hash_values["PL"].split(",").collect {|a| a.to_i}
    # use DP2 as other DP will overwrite this one
    hash_values["DP2"] = hash_values["DP"].to_i
    hash_values["GQ"] = hash_values["GQ"].to_f
    # now lets break up the info section too
    #infos = hash_values["INFO"].split(";").collect {|i| i.split("=")}
    #infos = infos.collect {|i| i.size == 1 ? i << true : i; i}
    #puts "#{infos.inspect}"
    #hash_infos = Hash[*infos.flatten]
    #puts "#{hash_infos.inspect}"
    #TODO: convert to floats / integers before adding to hash_values
    #hash_values.merge! hash_infos
    #puts "#{hash_values.inspect}"
    hash_values
  end

  def each
    header = []
    @file.each_line do |line|
      if line =~ /^##/
        next
      elsif line =~ /^#/
        header = line.chomp.gsub(/#/,"").split("\t")
      else
        raise "ERROR: header line not found" if header.empty?
        yield to_h(header, line)
      end
    end
  end

  def copy_if output_filename
    outfile = File.new(output_filename, 'w')
    header = []
    total_count = 0
    keep_count = 0
    kill_count = 0
    @file.each_line do |line|
      if line =~ /^##/
        outfile << line
        next
      elsif line =~ /^#/
        header = line.chomp.gsub(/#/,"").split("\t")
        outfile << line
        next
      else
        if header.empty?
          raise "ERROR: header line not found"
        end
        total_count += 1
        result = yield to_h(header, line)
        if result
          keep_count += 1
          outfile << line
        else
          kill_count += 1
        end
      end
    end
    outfile.close
    puts "kept   : #{keep_count} / #{total_count} (#{keep_count.to_f/total_count.to_f*100.0}%)"
    puts "removed: #{kill_count} / #{total_count} (#{kill_count.to_f/total_count.to_f*100.0}%)" 
  end

  def close
    @file.close
  end
end


