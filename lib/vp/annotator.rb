
class Annotator
  @@genome_databases = {"FruitFly" => "drosophila_melanogaster_core_57_513b"}

  attr_accessor :genome, :script_path

  def initialize(options = {})
    self.genome = options[:annotate]
    raise "ERROR: not valid genome: #{self.genome}. " unless valid_genome? self.genome
    self.script_path = find_script
    raise "ERROR: annotation script not found at: #{self.script_path}." unless self.script_path
  end

  def valid_genome? genome
    genome && @@genome_databases.keys.include?(genome)
  end

  def find_script
    annotate_script = File.join(File.dirname(__FILE__), "..","annotate", "annotateStrainSNVDiffs.pl")
    annotate_script = File.exists? annotate_script ? annotate_script : nil
    annotate_script
  end

  def database_for genome
    @@genome_databases[genome]
  end

  def annotate vcf_filename
    raise "ERROR: vcf file not found at: #{vcf_filename}" unless File.exists? vcf_filename
    database = database_for genome
    port = "5306"
    site = "ensembldb.ensembl.org"

    csv_filename = vcf_filename.split(".")[0..-2].join(".") + ".annotate.csv"
    log_filename = csv_filename + ".log"
    puts "Starting annotation on #{vcf_filename}"
    base_command = "#{self.script_path} "
    database_options = " #{@genome} #{database} #{port} #{site}"
    command = base_command + vcf_filename + database_options + " 1> #{csv_filename} 2> #{log_filename}"
    puts command
    system(command)
    csv_filename
  end

end
