#!/usr/bin/env ruby

require 'rubygems'
require 'bundler/setup'

require 'benchmark'
require 'logger'
include Benchmark

if ARGV.size != 2
  puts "ERROR invalid number of input parameters"
  puts "variant_annotator.rb <variant.file.vcf> <output.file.csv>"
  exit
end

vcfile = ARGV[0]
outfile = ARGV[1]

require 'ensembl'
require 'bio'

ActiveRecord::Base.logger = Logger.new(STDOUT)

include Ensembl::Core

puts "connecting to database"
DBConnection.connect("drosophila_melanogaster", 62)

class String
  def first_diff other_string
    diff_pos = 0
    self.each_char {|char| char != other_string[diff_pos] ? break : diff_pos += 1}
    diff_pos
  end
end

module Ensembl
  module Core
    
    # class ExonTranscript < DBConnection
    #    belongs_to :exon
    #    belongs_to :transcript
    # end
    
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
      
    end

    class Gene < DBConnection
      # same as above
      belongs_to :name_xref, :class_name => "Xref", :foreign_key => 'display_xref_id'
      #same as above
      def fast_display_name
        self.name_xref.display_label
      end
    end
    
    class Slice
      def fast_transcripts
        conditions = "seq_region_id = #{self.seq_region.id.to_s}"
        conditions += " AND seq_region_start >= #{self.start.to_s}"
        conditions += " AND seq_region_end <= #{self.stop.to_s}"
        Transcript.find(:all, :conditions => conditions, 
          :include => [:gene, :translation, {:exon_transcripts => :exon}, :transcript_stable_id, :seq_region])
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
    self.reference = reference.upcase
    self.length = self.reference.size
    self.stop = (self.start + self.length - 1).to_i
    self.quality = quality
    alleles = alleles.upcase
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
    @transcripts ||= self.slice.fast_transcripts
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

class Location
  attr_accessor :type, :strand, :position
  
  def initialize(type = :unknown, strand = 0, position = 0)
    self.type = type
    self.strand = strand
    self.position = position
  end
  
  def to_s
    location_string = "("
    # add one to position for display purposes
    location_string += "#{position + 1}:" unless type == :intron
    location_string += "#{strand}"
    location_string +=")" 
    location_string   
  end
end

class Impact
  attr_accessor :main, :detail, :message
  def initialize( main = nil, detail = nil, message = nil )
    self.main = main
    self.detail = detail
    self.message = message
  end
  
  def to_s
    result = ""
    result += "#{main}" if main
    result +=" (#{detail})" if detail
    result += " - #{message}" if message
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
      location = @transcript_location[transcript.id]
      location_info = location ? location.to_s : ""
      "#{transcript.fast_display_name}: #{location_info}"
    end
    transcript_strings.empty? ? "none" : transcript_strings.join("; ")
  end
  
  def impact_string
    impact_strings = @transcripts.map do |transcript| 
      impacts = @transcript_impacts[transcript.id]
      impact_info = impacts.inject("") {|info, impact| info + impact.to_s }
      "#{transcript.fast_display_name}: #{impact_info}"
    end
    impact_strings.empty? ? "none" : impact_strings.join("; ")
  end

  def to_s
    "#{variant_string}, #{gene_string}, #{transcript_string}, #{impact_string}\n"
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
      location = Location.new(:unkown,0)
      
      if biotype == 'protein_coding'
        gene = transcript.gene
        record.add_gene gene

        exclaim! "determining location of #{transcript.fast_display_name}"
        location = determine_location(variant, transcript)
        exclaim! "location done"
        
        exclaim! "determining impact of #{transcript.fast_display_name}"
        impacts = determine_impact(variant, transcript, location)
        exclaim! "impact done"
        
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
  
  def determine_location variant, transcript  
    # three possibilites for the location of the variant:
    #  * in the CDS
    #   ** exon is found and variant is within the start / stop of 
    #      the coding region
    #  * in an UTR
    #   ** exon is found but variant is outside of coding region
    #  * in a Intron
    #   ** no exon found    
    location = Location.new(:unknown, transcript.strand)
    
    variant_exon = transcript.exon_for_genomic_position(variant.start)
    coding_start = transcript.coding_region_genomic_start
    coding_stop = transcript.coding_region_genomic_end
    
    if variant_exon and ((variant.start > coding_start) and (variant.start < coding_stop))
      # in CDS
      exclaim! "variant in cds"
      location.type = :cds
      location.position = locate_mutation_in_cds(variant, transcript)
    elsif variant_exon and ((variant.start < coding_start) or (variant.start > coding_stop))
      # in UTR
      exclaim! "variant in UTR"
      location.type = :utr
      location.position = locate_mutation_in_cdna(variant, transcript)
    elsif !variant_exon
      #in intron
      exclaim! "variant in intron"
      location.type = :intron
    else
      puts "ERROR: invalid variant location"
      location.type = :unknown
    end
    location
  end
  
  def determine_impact variant, transcript, location
    impacts = Array.new
    begin
      case location.type
      when :cds
        impacts << analyze_in_cds(variant, transcript, location)
      when :utr
        impacts << analyze_in_utr(variant, transcript, location)
      when :intron
        impacts << analyze_in_intron(variant, transcript, location)
      else
        impacts << Impact.new("unknown")
      end
    rescue
      puts "ERROR analyzing #{transcript.fast_display_name} in #{location.type}"
      impacts << Impact.new("error")
    end
    impacts.flatten
  end
  
  def analyze_in_cds variant, transcript, location
    impacts = Array.new
    original_protein_sequence = transcript.protein_seq
    
    variant.mutant_alleles.each do |allele|      
      mutated_protein_sequence = mutate_transcript(allele, location.position, transcript)
      impacts << describe_protein_mutation(original_protein_sequence, mutated_protein_sequence)
    end
    impacts
  end

  def describe_protein_mutation original_sequence, mutated_sequence
    impact = Impact.new
    if !mutated_sequence.end_with? "*"
      # Stop codon has been lost. details for this impact 
      #  are the length of the translated sequence and what
      #  the stop codon has been changed to
      impact.main = "LS"
      impact.detail = "#{mutated_sequence.size}, #{original_sequence[-1..-1]}->#{mutated_sequence[-1..-1]}"
    elsif original_sequence == mutated_sequence
      # Synonymous resulting mutated protein
      # we do not provide more info
      impact.main = "synonymous"
    elsif ((stop_codon_position = mutated_sequence =~ /\*/) != (mutated_sequence.size - 1))
      impact.main = "nonsense"
      impact.detail = "#{stop_codon_position}:#{mutated_sequence.size}"
    else
      impact.main = "other"
      # find where they don't match
      if original_sequence.size != mutated_sequence.size
        puts "ERROR: original and mutated protein length do not match" 
      else
        diff_pos = original_sequence.first_diff mutated_sequence
        old_peptide = original_sequence[diff_pos]
        new_peptie = mutated_sequence[diff_pos]
        diff_pos += 1 # for display
        impact.detail = "#{old_peptide}:#{new_peptie} #{diff_pos}"
      end
    end
    impact
  end
  
  def locate_mutation_in_cds variant, transcript
    transcript.genomic2cds(variant.start)
  end
  
  def locate_mutation_in_cdna variant, transcript
    transcript.genomic2cdna(variant.start)
  end

  def mutate_transcript raw_mutation, mutation_location, transcript
    # cds_seq returns the coding sequence of the transcript
    #  the concatenated sequence of all exons minus the UTRs
    sequence = transcript.cds_seq
    
    # flip the mutation if the strand is reversed
    # TODO: not sure if this is necessary. Was in the old code
    mutation = (transcript.strand == 1) ? raw_mutation : 
      raw_mutation.gsub(/(a|t|c|g)/, {'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c'}) 
      
    mutated_sequence = sequence
    mutated_sequence[mutation_location] = mutation

    Bio::Sequence::NA.new(mutated_sequence).translate
  end
  
  def analyze_in_utr variant, transcript, location
    impact = Impact.new "utr"
  end
  
  # mostly just copied from the old perl script
  def analyze_in_intron variant, transcript, location
    impact = Impact.new "intron"
    introns = transcript.introns
    introns.each do |intron|
      if variant.start == intron.seq_region_start or variant.stop == intron.seq_region_start
        exon = transcript.strand == 1 ? intron.previous_exon : intron.next_exon
        impact.message = "SJ;Exon: #{exon.seq_region_id}"
        break
      elsif variant.start == intron.seq_region_end or variant.stop == intron.seq_region_end
        exon = transcript.strand == 1 ? intron.next_exon : intron.previous_exon
        impact.message = "SJ;Exon: #{exon.seq_region_id}"
        break
      elsif variant.start > intron.seq_region_start and variant.stop < intron.seq_region_end
        next_exon = transcript.strand == 1 ? intron.next_exon : intron.previous_exon
        previous_exon = transcript.strand == 1 ? intron.previous_exon : intron.next_exon
        
        next_exon_distance = (variant.stop - intron.next_exon.seq_region_start).abs
        previous_exon_distance = (variant.start - intron.previous_exon.seq_region_end).abs

        message = "E: #{previous_exon.seq_region_id} #{previous_exon_distance} bp"
        message += "E: #{next_exon.seq_region_id} #{next_exon_distance} bp"
        
        if next_exon_distance == 20
          # a branch point is 20bps away from the AG side of the intron (aka the end of the intron)
          message += " BP"
        end
        impact.message = message
        break
      end      
    end
    impact
  end

  def out
    exclaim! "writing out to file"
    header = "chromosome, variant_start, variant_stop, reference, mutation, quality, genes, transcripts(pos:strand), impact"
    records_out = @records.join()
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

