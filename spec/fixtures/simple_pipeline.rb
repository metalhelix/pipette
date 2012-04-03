
class SimplePipeline < Pipeline
  name "simple"
  description "a simple pipeline"
  step :step_1 do
    input :input_1
    output :output_1 do |input| "#{input[:input_1]}_out" end
    run do |inputs, outputs|
      "run_#{inputs[:input_1]}_#{outputs[:output_1]}"
    end
  end
end

class TwoStepPipeline < Pipeline
  default_steps :step_1, :step_2
  step :step_1 do
    input :input_1
    output :output_1 do |input| "#{input[:input_1]}_out" end
    run do |inputs, outputs|
      "step_1_#{inputs[:input_1]}"
    end
  end
  step :step_2 do
    input :output_1
    output :output_2 do |input| "#{input[:output_1]}_out" end
    run do |inputs, outputs|
      "step_2_#{inputs[:output_1]}"
    end
  end
end

require 'optparse'
class OptionsPipeline < Pipeline
  options do
    opts = OptionParser.new do |o|
      o.on('-i', '--input BAM_FILE', 'REQUIRED - Input BAM file to call SNPs on') {|b| options[:input] = b}
      o.on('-r', '--reference FA_FILE', 'REQUIRED - Reference Fasta file for genome') {|b| options[:reference] = b}
      o.on('-o', '--output PREFIX', 'Output prefix to use for generated files') {|b| options[:output] = b}
      o.on('-h', '--help', 'Displays help screen, then exits') {puts o; exit}
    end
    opts
  end
  step :test_input do

  end
end
