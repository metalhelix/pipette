class Step
  attr_reader :name, :inputs, :outputs
  attr_accessor :run_block

  def initialize name #:nodoc:
    @name = name
    @inputs = []
    @outputs = []
    @output_blocks = {}
    @run_block = nil
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

  # Provide the name of an input value that this step is dependent on
  # Inputs are passed to the run block when executing the step.
  # If inputs are not all provided, step and thus pipeline cannot and
  # will not run.
  # ==== Parameters
  # name<Symbol>:: Name of the input that step needs to run
  def input name
    @inputs << name.to_sym
  end

  # Provide the name of an output that this step will produce when
  # it is run. The block passed in is given the inputs to the step so
  # that the step can generate the return value for this output parameter
  # This allows the output content to still be determined by the step that is
  # producing the output but also allows the step to be skipped and still be
  # able to recover the output content it should produce.
  # ==== Parameters
  # name<Symbol>:: name of the output the step will generate
  # block<Block>:: block that will determine output content given input content
  def output name, &block
    name = name.to_sym
    @outputs << name
    @output_blocks[name] = block if block
  end

  # Provide the run block that will be executed when this step is run.
  # This block is the actual functionality of the step. It is passed
  # the inputs and outputs as parameters upon execution
  # ==== Parameters
  # run_block<Block>:: code to execute when running this step
  def run &run_block
    self.run_block = run_block
  end

  # Given inputs, execute the blocks provided to the output method to
  # generate ouput content for the step
  # ==== Parameters
  # inputs<Hash>:: inputs provided to the step
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
  # inputs<Hash>:: inputs hash with :input_name => input value
  # ==== Returns
  # Array[Symbol]
  def missing_inputs inputs
    self.inputs.select {|input| !inputs.include?(input)}
  end

  # Executes the run block provided by the 'run' method when creating the step
  # Raises exceptions if no run block was given to this step
  # or if there are missing inputs parameters in the inputs
  # ==== Parametesr
  # Inputs<Hash>:: inputs hash with :input_name => "input_value" format
  # ==== Returns
  # Result of calling the Run block
  def call_run_block(inputs, outputs)
    raise "No call section provided for #{name}" unless @run_block
    missing = missing_inputs(inputs)
    raise "Missing inputs: #{missing.join(",")}" unless missing.empty? 
    @run_block.call(inputs, outputs) if @run_block
  end
end
