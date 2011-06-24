
require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe Step do
  before(:each) do
    @step = Step.new :step
  end

  it "should have a name" do
    @step.name.should == :step
    s = Step.new "string name"
    s.name.should == "string name"
  end

  it "should have no inputs or outputs initially" do
    @step.inputs.empty?.should == true
    @step.outputs.empty?.should == true
  end

  it "should take input and output names" do
    @step.input :input_1
    @step.input :input_2
    @step.inputs.size.should == 2
    @step.output :output_1
    @step.outputs.size.should == 1
    @step.output :output_2
    @step.outputs.size.should == 2
    @step.output :output_3 do |input| "output_3" end
    @step.outputs.size.should == 3
  end

  it "should evaluate output parameters given inputs" do
    inputs = {:input_1 => "input_1"}
    @step.input :input_1
    @step.output :output_1 do |input| "#{input[:input_1]}.out" end
    @step.outputs_given(inputs).size.should == 1
    @step.outputs_given(inputs)[:output_1].should == "input_1.out"
  end

  it "should run its run block" do
    @step.run_block = lambda{ |input,output| "run_#{input}" }
    @step.run_block.should_not == nil
    output = @step.call_run_block("input", "output").should == "run_input"
  end

  it "should have inputs and outputs populated in run" do
    input = {}
    input[:input_1] = "111"
    input[:input_2] = "222"
    @step.input :input_1
    @step.input :input_2
    @step.output :output_1 do |input| "#{input[:input_1]}_output" end
    output = @step.outputs_given(input)
    @step.run_block = lambda{ |input,output| "run_#{input[:input_1]}_#{input[:input_2]}_#{output[:output_1]}" }
    output = @step.call_run_block(input,output)
    output.should == "run_111_222_111_output"
  end

end
