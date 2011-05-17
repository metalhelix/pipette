
class Annotator
  attr_accessor :database, :jar_path, :config_path

  def initialize(options = {})
    self.jar_path = options[:snpeff]
    raise "ERROR: snpEff jar not found at: #{self.jar_path}." unless File.exists?(self.jar_path)
    self.config_path = options[:snpeff_config]
    raise "ERROR: snpEff config not found at: #{self.config_path}." unless File.exists?(self.config_path)
    self.database = options[:annotate]
    raise "ERROR: not valid annotation database: #{self.database}. " unless valid_database? self.database
  end

  def valid_database? database
    true
  end

  def annotate vcf_filename
    raise "ERROR: vcf file not found at: #{vcf_filename}" unless File.exists? vcf_filename

    base_name = vcf_filename.split(".")[0..-2].join(".")
    output_filename = base_name + ".snpeff.txt"
    stats_filename = base_name + ".snpeff.summary.html"
    puts "Starting annotation on #{vcf_filename}"
    command = "java -Xmx4g -jar #{self.jar_path} -c #{config_path}"
    command += " -no-downstream -no-upstream -ud 0 -vcf4"
    command += " -stats #{stats_filename}"
    command += " #{self.database} #{vcf_filename} > #{output_filename}"

    puts command
    system(command)
    output_filename
  end
end
