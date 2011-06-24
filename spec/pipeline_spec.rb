require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe Pipeline do
  describe "creation" do
    before(:each) do
      @p = Pipeline.new
    end
    it "should initialize" do
      @p.steps.should == []
    end
    it "should have defaults steps" do
      @p.default_steps.should == nil
      steps = [:a, :b]
      @p.default_steps = steps
      @p.default_steps.should == steps
    end
  end
  describe "with steps" do
    before(:each) do
      @p = Pipeline.new
      s1 = Step.new :step_1
      s1.input :input_1
      s1.output :output_1
      s1.run_block = lambda { "run" }
      @p.steps << s1
      s2 = Step.new :step_2
      s2.input :output_1
      s2.input :missing_1
      @p.steps << s2
    end
    it "should have two steps" do
      @p.steps.size.should == 2
    end
    it "should find missing input parameters" do
      input = {:not_input => "not"}
      result = @p.missing_step_inputs input
      result.should == [[:step_1, :input_1],[:step_2, :missing_1]]
      content = capture(:stdout){ @p.run(input) }
      content.should =~ /step_1.*missing.*input_1/
      content.should =~ /step_2.*missing.*missing_1/
      input = {:input_1 => "input"}
      result = @p.missing_step_inputs input
      result.should == [[:step_2, :missing_1]]
      @p.steps[0].output :missing_1
      result = @p.missing_step_inputs input
      result.should == []
    end
  end
  describe "simple pipeline" do
    before(:each) do
      @simple = SimplePipeline.pipeline
    end
    it "should have one step" do
      @simple.steps.size.should == 1
    end
    it "should have input and output in its step" do
      @simple.steps[0].inputs.size.should == 1
      @simple.steps[0].inputs.include?(:input_1).should == true
      @simple.steps[0].outputs.size.should == 1
      @simple.steps[0].outputs.include?(:output_1).should == true
    end
    it "should run its run section" do
      inputs = {:input_1 => 'input'}
      result = @simple.run inputs
      result.should == ["run_input_input_out"]
    end
  end
  describe "two step pipeline" do
    before(:each) do
      @two_step = TwoStepPipeline.pipeline
    end
    it "should have two steps" do
      @two_step.steps.size.should == 2
    end
    it "should have two default steps" do
      @two_step.default_steps.size.should == 2
    end
  end
end
