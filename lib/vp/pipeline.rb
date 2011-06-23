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

  def input value
    @inputs << value.to_sym
  end

  def output value, &block
    value = value.to_sym
    @outputs << value
    @output_blocks[value] = block if block
  end

  def outputs_given inputs
    @output_values = {}
    @outputs.each do |output|
      @output_values[output] = @output_blocks[output].call(inputs) if @output_blocks[output]
    end
    @output_values
  end

  def run(inputs, outputs)
    @run_block.call(inputs, outputs) if @run_block
  end
end

class Pipeline

  def self.step name, &block
    @current_step = Step.new(name)
    block.call
  end

  def self.input name
    @current_step.input name
  end

  def self.output name, &block
    @current_step.output(name, &block)
  end

  def self.run &run_block
    @current_step.run_block = run_block
    @@steps ||= []
    @@steps << @current_step
    @current_step = nil
  end

  def initialize
    @steps = []
    @current_step = nil
  end

  def run options
    check_steps options[:steps]

    options[:steps].each do |step|
      results = self.send(step.to_sym, options)
      options.merge! results
    end
  end

  def execute command
    report command
    result = system(command)
  end

  def report status
    puts "#{Time.now} - " + status 
  end

  def check_steps steps
    raise "ERROR: no steps given" unless steps
    raise "ERROR: steps not array" unless steps.respond_to? :each
    steps.each do |step|
      if !self.respond_to?(step.to_sym)
        raise "ERROR: step #{step} not valid"
      end
    end
  end
end
