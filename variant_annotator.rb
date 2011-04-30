#!/usr/bin/env ruby

require 'benchmark'
include Benchmark

if ARGV.size != 2
  puts "ERROR invalid number of input parameters"
  puts "variant_annotator.rb <variant.file.vcf> <output.file.csv>"
  exit
end

vcfile = ARGV[0]
outfile = ARGV[1]

puts "require ensembl gem"
require 'ensembl'
include Ensembl::Core

puts "require bio gem"
require 'bio'

puts "connecting to database"
DBConnection.connect("drosophila_melanogaster", 62)

module Ensembl
  module Core
    class Transcript < DBConnection
      # This directly associates the Xref for the transcripts
      # display label. This way, the display_label for the
      # transcript is pulled from the database when we query
      # for the transcript
      belongs_to :name_xref, :class_name => "Xref", :foreign_key => 'display_xref_id'

      # Additional method to use our name_xref association to
      # get the name, instead of the slow way the ensembl gem
      # is doing it now
      def fast_display_name
        self.name_xref.display_label
      end
      #has_many :exon_transcripts, :include => :exon
      #has_many :exons, :through => :exon_transcripts, :order => "rank ASC"
    end

    class Gene < DBConnection
      # same as above
      belongs_to :name_xref, :class_name => "Xref", :foreign_key => 'display_xref_id'

      #same as above
      def fast_display_name
        self.name_xref.display_label
      end
    end
  end
end

# Represents a variant from a line in the vcf file.
# holds all the data pulled from that line and allows us to 
# work with this info in a object oriented way
class Variant
  attr_accessor :start, :stop, :chromosome, :reference, :alleles, :length, :quality, :mutant_alleles

  def initialize chromosome, start, reference, quality, alleles
    self.chromosome = chromosome
    self.start = start.to_i
    self.reference = reference.downcase
    self.length = self.reference.size
    self.stop = (self.start + self.length - 1).to_i
    self.quality = quality
    alleles = alleles.downcase
    if alleles.kind_of? String
      self.alleles = alleles.split(",")
    else
      self.alleles = alleles
    end
    self.mutant_alleles = self.alleles.reject {|allele| allele == self.reference}
  end

  def to_s
    "#{chromosome}, #{start}, #{stop}, #{reference}, #{alleles.join(";")}, #{quality}"
  end

  def self.from_vcf vcf_line
    chr, start_pos, dummy, ref_chr, mut_allels, quality = vcf_line.split("\t")
    Variant.new(chr, start_pos, ref_chr, quality, mut_allels)
  end

end

# Chromosome class is used to 
# store aggregate data from the database
# this allows us to not hit the database for
# every new variant we are annotating
class Chromosome
  attr_accessor :name
  
  def initialize(name)
    self.name = name
    @transcripts = nil
    @slice = nil
  end

  def transcripts_overlap(variant)
    self.transcripts.select {|tran| tran.start <= variant.start and variant.stop <= tran.stop}
  end

  def slice
    # This will lazy load the slice representing the 
    # entire chromosome, and store it in the member variable
    # so that next time we will just use the data in @slice
    # and not have to hit the database again
    @slice ||= Slice.fetch_by_region('chromosome', self.name)
  end

  def transcripts
    # Here we lazy load all the transcripts for the
    # chromosome and cache them in @transcripts so we
    # only hit the database once
    @transcripts ||= self.slice.transcripts(:include => [ :exon_transcripts, :seq_region ] )
  end
end

class Genome
  @@chromosomes = {}

  def self.chromosome name
    # if we have acquired this chromosome
    # before, then get it out of the hash
    # else, create a new one with this name now
    @@chromosomes[name] ||= Chromosome.new(name)

    @@chromosomes[name]
  end
end

class Impact
  attr_accessor :main, :detail, :message
  def initialize( main = nil, detail = nil, message = nil)
    self.main = main
    self.detail = detail
    self.message = message
  end
  
  def to_s
    result = ""
    result += "#{main}" if main
    result +=" (#{detail})" if detail
    result += " #{message}" if message
    result
  end
end

# Maintains data that will be used to 
# annotate the variant. knows how to 
# output the data in the way we want it
class Record
  def initialize variant
    @variant = variant
    @genes = Array.new 
    @transcripts = Array.new
    @transcript_impacts = Hash.new
    @transcript_location = Hash.new
  end

  def add_gene gene
    puts "ERROR: not a gene: #{gene}" unless gene.kind_of? Ensembl::Core::Gene
    @genes << gene unless @genes.include? gene
  end
  
  def add_transcript transcript, location, impacts
    puts "ERROR: not a transcript: #{transcript}" unless transcript.kind_of? Ensembl::Core::Transcript
    unless @transcripts.include? transcript
      @transcripts << transcript 
      @transcript_location[transcript.id] = location
      @transcript_impacts[transcript.id] = impacts unless impacts.empty?
    end
  end

  def variant_string
    @variant.to_s
  end

  def gene_string
    @genes.empty? ? "none" : @genes.map {|gene| "#{gene.fast_display_name}"}.join("; ")
  end

  def transcript_string
    transcript_strings = @transcripts.map do |transcript| 
      impacts = @transcript_impacts[transcript.id]
      additional_info = impacts.inject("") {|info, impact| info + impact.to_s }
      "#{transcript.fast_display_name}:#{additional_info}(#{transcript.seq_region_start}:#{transcript.seq_region_strand})"
    end
    transcript_strings.empty? ? "none" : transcript_strings.join("; ")
  end

  def to_s
    "#{variant_string}, #{gene_string}, #{transcript_string}\n"
  end
end


# performs the actual annotation of each
# variant found in the vcf file
class Annotator
  def initialize input_file, output_file_name
    @input_file = input_file
    @output_file_name = output_file_name
    @records = []
    @verbose = true
  end

  def exclaim! str
    puts str if @verbose
  end

  def run 
    file = File.open(@input_file, 'r')
    file.each_line do |line|
      # skip comment lines
      next if line =~ /^#/
      variant = Variant.from_vcf line
      # skip if there aren't really any 
      # variant alleles in the variant
      next if variant.mutant_alleles.empty?
      # valid variant - lets annotate it!
      annotate_variant variant
    end
    file.close
  end

  def annotate_variant variant
    # record will store relevant info about
    # variant we are annotating
    record = Record.new variant
    # acquire the chromosome this variant is on 
    # from our Genonme 'database' 
    # This allows us to load all transcripts 
    # on a chromosome at one time
    chromosome = Genome.chromosome(variant.chromosome)

    # acquire transcripts that this variant overlaps
    exclaim! "getting overlapping transcripts"
    transcripts = chromosome.transcripts_overlap variant 
    exclaim! "num transcripts: #{transcripts.size}"

    # loop through each transcript to gather more info
    transcripts.each do |transcript|

      biotype = transcript.biotype
      
      impacts = Array.new
      location = :unknown
      
      if biotype == 'protein_coding'
        gene = transcript.gene
        record.add_gene gene

        location = determine_location(variant, transcript)
        impacts = determine_impact(variant, transcript, location)
        
      else
        impacts << Impact.new(biotype)
      end
      
      exclaim! "recording transcript"
      record.add_transcript transcript, location, impacts
    end

    exclaim! "vcf line done. adding to records"
    @records  << record.to_s
    exclaim! "done"
  end
  
  def determine_impact variant, transcript, location
    impacts = Array.new
    case location
    when :exon
      impacts << analyze_in_exon(variant, transcript)
    when :utr
      impacts << analyze_in_utr(variant, transcript)
    when :intron
      impacts << analyze_in_intron(variant, transcript)
    else
      impacts << Impact.new("unknown")
    end
    impacts.flatten
  end

  def determine_location variant, transcript  
    # three possibilites for the location of the variant:
    #  * in the CDS
    #   ** exon is found and variant is within the start / stop of 
    #      the coding region
    #  * in an UTR
    #   ** exon is found but variant is outside of coding region
    #  * in a Intron
    #   ** no exon found
    location = :unknown
    variant_exon = transcript.exon_for_genomic_position(variant.start)
    coding_start = transcript.coding_region_genomic_start
    coding_stop = transcript.coding_region_genomic_end
    
    if variant_exon and ((variant.start > coding_start) and (variant.start < coding_stop))
      # in CDS
      exclaim! "variant in exon"
      location = :exon
    elsif variant_exon and ((variant.start < coding_start) or (variant.start > coding_stop))
      # in UTR
      exclaim! "variant in UTR"
      location = :utr
    elsif !variant_exon
      #in intron
      exclaim! "variant in intron"
      location = :intron
    else
      puts "ERROR: invalid variant location"
      location = :unknown
    end
    location
  end
  
  def analyze_in_exon variant, transcript
    impacts = Array.new
    original_sequence = transcript.protein_seq
    variant.mutant_alleles.each do |allele|
      mutated_sequence = mutate_transcript(allele, variant.start, transcript)
      impacts << analyze_exon_mutation(original_sequence, mutated_sequence)
    end
    impacts
  end

  def analyze_exon_mutation original_sequence, mutated_sequence
    impact = Impact.new
    if !mutated_sequence.end_with? "*"
      impact.main = "lost stop codon"
      impact.detail = "1,1"
    elsif original_sequence == mutated_sequence
      impact.main = "synonymous"
    elsif (mutated_sequence =~ /\*/) != (mutated_sequence.size - 1)
      impact.main = "nonsense"
    else
      impact.main = "other"
    end
    impact
  end

  def mutate_transcript raw_mutation, raw_start, transcript
    # cds_seq returns the coding sequence of the transcript
    #  the concatenated sequence of all exons minus the UTRs
    sequence = transcript.cds_seq
    mutation_location = transcript.genomic2cds(raw_start)
    # flip the mutation if the strand is reversed
    # TODO: not sure if this is necessary. Was in the old code
    mutation = (transcript.strand == 1) ? raw_mutation : 
      raw_mutation.gsub(/(a|t|c|g)/, {'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c'}) 
      
    mutated_sequence = sequence
    mutated_sequence[mutation_location] = mutation

    Bio::Sequence::NA.new(mutated_sequence).translate
  end
  
  def analyze_in_utr(variant, transcript)
    impact = Impact.new
    
  end
  
  def analyze_in_intron(variant, transcript)
    impact = Impact.new
  end

  def out
    exclaim! "writing out to file"
    header = "chromosome, variant_start, variant_stop, reference, mutation, quality, genes, transcripts(pos:strand), location, impact, go terms"
    records_out = @records.join("\n")
    output_file = File.new(@output_file_name, 'w')
    output_file << header << "\n"
    output_file << records_out
    output_file.close
  end
end

time = Benchmark.measure do
annotator = Annotator.new vcfile, outfile
annotator.run
annotator.out
end

puts time

#Benchmark.bm(7) do |x|
#  x.report("run:") { annotator.run }
#  x.report("out:") { annotator.out }
#end

