class Step
  attr_reader :name, :inputs, :outputs
  attr_accessor :run_block

  def initialize name
    @name = name
    @inputs = []
    @outputs = []
    @output_blocks = {}
    @run_block = nil
  end

  def evaluate &block
    # trick from
    # http://www.dan-manges.com/blog/ruby-dsls-instance-eval-with-delegation
    @self_before_instance_eval = eval "self", block.binding
    instance_eval(&block)
  end

  def method_missing(method, *args, &block)
    # trick from
    # http://www.dan-manges.com/blog/ruby-dsls-instance-eval-with-delegation
    @self_before_instance_eval.send method, *args, &block
  end

  def input value
    @inputs << value.to_sym
  end

  def output value, &block
    value = value.to_sym
    @outputs << value
    @output_blocks[value] = block if block
  end

  def run &run_block
    self.run_block = run_block
  end

  def outputs_given inputs
    @output_values = {}
    @outputs.each do |output|
      @output_values[output] = @output_blocks[output].call(inputs) if @output_blocks[output]
    end
    @output_values
  end

  # Returns an array of inputs that the step requires but are NOT present in the
  # 'inputs' input parameter.
  # ==== Parameters
  # inputs<Hash>:: inputs hash with input_name => input value
  # ==== Returns
  # Array[Symbol]
  def missing_inputs inputs
    self.inputs.select {|input| !inputs.include?(input)}
  end

  # Executes the run block provided by the 'run' method when creating the step
  # Raises exceptions if no run block was given to this step
  # or if there are missing inputs parameters in the inputs
  def call_run_block(inputs, outputs)
    raise "No call section provided for #{name}" unless @run_block
    missing = missing_inputs(inputs)
    raise "Missing inputs: #{missing.join(",")}" unless missing.empty? 
    @run_block.call(inputs, outputs) if @run_block
  end
end

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

  def run inputs
    check_steps inputs
    outputs = {}
    results = []
    @steps.each do |step|
      step_output = step.outputs_given inputs
      outputs.merge! step_output
      results << step.call_run_block(inputs,outputs)
    end
    results
  end

  def check_steps inputs

  end

  def execute command
    report command
    result = system(command)
  end

  def report status
    puts "#{Time.now} - " + status 
  end
end
