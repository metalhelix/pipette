require 'vp/step'

class Pipeline
  attr_accessor :steps, :default_steps

  def self.step name, &block
    current_step = Step.new(name)
    current_step.evaluate(&block)
    pipeline.steps << current_step
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
    @steps = []
    @default_steps = nil
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
  def run inputs
    run_step_names = run_steps inputs
    puts "running steps: #{run_step_names.join(", ")}"
    missing = missing_step_inputs inputs
    output_missing(missing) and return unless missing.empty?
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
    puts "all steps    : #{all_step_names.join(", ")}"
    run_step_names = default_steps
    if inputs and inputs[:steps]
      input_steps = [inputs[:steps]].flatten.collect {|s| s.strip.downcase.to_sym}
      puts "input steps options: #{input_steps.join(', ')}"
      run_step_names = all_step_names.select {|s| input_steps.include? s}
    end

    if run_step_names.nil?
      puts "defaulting to all steps"
      run_step_names = all_step_names
    end

    run_step_names
  end

  # Prints the provided list of missing input
  # parameters to standard out.
  # ==== Parameters
  # missing<Array>:: Array of arrays of missing parameters
  # for each step
  def output_missing(missing)
    puts "ERROR: Missing input parameters for steps"
    output = missing.collect {|step, miss| "#{step} is missing #{miss}"}
    output.each {|o| puts o}
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
