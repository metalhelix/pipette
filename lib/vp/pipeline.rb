require 'vp/step'

class Pipeline
  attr_accessor :steps, :default_steps

  def self.step name, &block
    current_step = Step.new(name)
    current_step.evaluate(&block)

    pipeline.steps << current_step
  end

  def self.pipeline
    @pipeline ||= self.new
    @pipeline
  end

  def self.default_steps *steps
    pipeline.default_steps = steps
  end

  def self.inherited(subclass)
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
    missing = missing_step_inputs inputs
    output_missing(missing) and return unless missing.empty?
    results = []
    @steps.each do |step|
      step_output = step.outputs_given inputs
      results << step.call_run_block(inputs,step_output)
      inputs.merge! step_output
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
