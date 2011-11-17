
class SnpEff
  def initialize
  end

  #TODO: check if database is actually valid
  def valid_database? database
    true
  end

  def check_options options
    options[:snpeff] = File.expand_path(options[:snpeff])
    raise "ERROR: snpEff jar not found at: #{options[:snpeff]}." unless File.exists?(options[:snpeff])
    options[:snpeff_config] = File.expand_path(options[:snpeff_config])
    raise "ERROR: snpEff config not found at: #{options[:snpeff_config]}." unless File.exists?(options[:snpeff_config])
    raise "ERROR: not valid annotation database: #{options[:annotate]}. " unless valid_database? options[:annotate]
  end

  def run vcf_filename, options
    check_options options
    jar_path = options[:snpeff]
    database = options[:annotate]
    config_path = options[:snpeff_config]
    raise "ERROR: vcf file not found at: #{vcf_filename}" unless File.exists? vcf_filename

    base_name = vcf_filename.split(".")[0..-2].join(".")
    output_filename = base_name + ".snpeff.txt"
    stats_filename = base_name + ".snpeff.summary.html"
    puts "Starting annotation on #{vcf_filename}"
    command = "java -Xmx4g -jar #{jar_path} -c #{config_path}"
    command += " -no-downstream -no-upstream -ud 0"
    command += " -stats #{stats_filename}"
    command += " #{database} #{vcf_filename} > #{output_filename}"

    puts command if options[:verbose]
    system(command)
    results = {:snpeff_output => output_filename, :snpeff_summary => stats_filename}
    results
  end
end
