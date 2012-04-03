require 'optparse'
require 'yaml'
require 'fileutils'

require 'pipette/step'

class Pipeline
  attr_accessor :name, :description, :steps, :default_steps, :options, :options_parser

  # DSL method for defining a step
  # takes the name of the step and the
  # block to evaluate to become a step
  def self.step name, &block
    current_step = Step.new(name)
    current_step.evaluate(&block)
    pipeline.steps << current_step
  end

  # DSL method for giving a pipeline a
  # name
  def self.name name
    # replace spaces with underscores
    pipeline.name = name.gsub(/\s+/,"_")
  end

  # DSL method for giving a pipeline a
  # description
  def self.description description
    pipeline.description = description
  end

  # DSL method for defining an options parser
  # takes a block which should return an OptionsParser
  # instance to parse with the input arguments
  def self.options &block
    pipeline.options_parser = pipeline.evaluate(&block)
  end

  # Returns instance of pipeline.
  # Should be used to get the pipeline
  # created from pipeline file (instead of new)
  def self.pipeline
    @pipeline ||= self.new
    @pipeline
  end

  # provide pipeline with a list of
  # step names. These will be the
  # steps run if specific steps to run
  # are not provided on the command line
  # A value of nil (the default value
  # for default_steps indicates that
  # all steps of the pipeline will be
  # run by default
  def self.default_steps *steps
    pipeline.default_steps = steps
  end

  def self.inherited(subclass) #:nodoc:
    # could be cool to use
    # callback that is executed when a new
    # sub-class is created
  end

  def initialize
    @name = "default"
    @description = ""
    @steps = []
    @default_steps = nil
    @options = {}
    @options_parser = nil
  end

  def default_options
    default_options_file = File.expand_path(File.join(File.dirname(__FILE__), "..", "..", "config", "#{self.name}_config.yml"))
    if File.exists? default_options_file
      self.options = Hash[YAML::load(open(default_options_file)).map {|k,v| [k.to_sym, v]}]
    else
      puts "WARNING: no default configuration found in #{default_options_file}"
    end
  end

  def evaluate &block #:nodoc:
    # trick from
    # http://www.dan-manges.com/blog/ruby-dsls-instance-eval-with-delegation
    @self_before_instance_eval = eval "self", block.binding
    instance_eval(&block)
  end

  def method_missing(method, *args, &block) #:nodoc:
    # trick from
    # http://www.dan-manges.com/blog/ruby-dsls-instance-eval-with-delegation
    # see above 'evaluate' as well
    @self_before_instance_eval.send method, *args, &block
  end

  # runs the steps of the pipeline with the
  # inputs provided. Each step gets the
  # combination of the starting input and the
  # the outputs of the previous steps.
  # Steps can overwrite old inputs with new
  # outputs. run checks if required inputs
  # for each step is satisfied
  # ==== Parameters
  # inputs<Hash>:: Initial set of inputs provided
  # by the user.
  def run inputs = {}
    inputs.merge!(self.options)
    output_inputs(inputs)
    run_step_names = run_steps inputs
    puts "running steps: #{run_step_names.join(", ")}"
    missing = missing_step_inputs inputs
    unless missing.empty?
      output_missing(missing)
      return
    end
    results = []
    @steps.each do |step|
      puts "step: #{step.name}"
      step_output = step.outputs_given inputs
      if run_step_names.include? step.name
        results << step.call_run_block(inputs,step_output)
        inputs.merge! step_output
      else
        inputs.merge! step_output
      end
    end
    results
  end

  def parse_input input_args
    if self.options_parser
      self.options_parser.parse!(input_args)
    end
    options
  end

  # Returns an array of arrays, one for
  # each input missing for each step.
  # The output from the previous step is used
  # as input the next step, so these step dependent
  # inputs will not count towards the missing steps.
  # If there are no missing steps, an empty array is
  # returned
  # ==== Parameters
  # inputs<Hash>:: Initial set of inputs provided by user
  # ==== Returns
  # Array of arrays. Elements of the internal array have
  # the following content:
  # [step_name, missing_input]
  def missing_step_inputs inputs
    inputs_names = inputs.keys
    inputs_names ||= []
    results = []
    @steps.each do |step|
      output_names = step.outputs
      missing = step.missing_inputs inputs_names
      inputs_names += output_names
      missing.each do |miss|
        results << [step.name, miss]
      end
    end
    results
  end

  def run_steps inputs = nil
    all_step_names = @steps.collect {|s| s.name}
    # puts "all steps    : #{all_step_names.join(", ")}"
    run_step_names = default_steps
    if inputs and inputs[:steps]
      input_steps = [inputs[:steps]].flatten.collect {|s| s.strip.downcase.to_sym}
      # puts "input steps options: #{input_steps.join(', ')}"
      run_step_names = all_step_names.select {|s| input_steps.include? s}
    end

    if run_step_names.nil?
      puts "defaulting to all steps"
      run_step_names = all_step_names
    end

    run_step_names
  end

  def output_inputs(inputs)
    puts "inputs used:"
    inputs.each do |input, value|
      puts "  #{input} => #{value}"
    end
    puts ""
  end

  # Prints the provided list of missing input
  # parameters to standard out.
  # ==== Parameters
  # missing<Array>:: Array of arrays of missing parameters
  # for each step
  def output_missing(missing)
    puts ""
    puts "ERROR: Missing input parameters for steps:"
    output = missing.collect {|step, miss| "  #{step} is missing #{miss}"}
    output.each {|o| puts o}
    puts ""
  end

  # Helper function that executes a command line command as well as
  # reporting this command to output
  def self.execute command
    report command
    result = system(command)
  end

  # Output status parameter as well as other logging information
  def self.report status
    puts "#{Time.now} - " + status
  end
end
