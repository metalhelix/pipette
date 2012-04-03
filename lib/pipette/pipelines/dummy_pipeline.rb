require 'pipette/pipeline'

class DummyPipeline < Pipeline
  name "dummy"
  description "simply prints out inputs and outputs"
  options do
    opts = OptionParser.new do |o|
      o.on('-i', '--input BAM_FILE', 'REQUIRED - Input BAM file to call SNPs on') {|b| options[:input] = b}
      o.on('-r', '--reference FA_FILE', 'REQUIRED - Reference Fasta file for genome') {|b| options[:reference] = b}
    end
    opts
  end

  step :print_input do
    run do |input, output|
      puts "input:"
      puts input.inspect
      puts "output:"
      puts output.inspect
    end
  end
end
