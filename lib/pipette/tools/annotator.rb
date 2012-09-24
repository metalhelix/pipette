require 'pipette/tools/snpeff'
require 'pipette/tools/snpeff_summarizer'

class Annotator
  def run input_file, options
    snpeff = SnpEff.new
    results = snpeff.run input_file, options
    summer = SnpeffSummarizer.new
    options[:raw_vcf_file] = input_file
    puts "summarize: #{results[:snpeff_output]}"
    summary_file = summer.run results[:snpeff_output], options
    summary_file
  end
end

