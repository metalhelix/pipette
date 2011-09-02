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

    it "should have default running steps by default" do
      input = {:input_1 => "in"}
      @p.run_steps(input).size.should == 2
      @p.default_steps = [:step_1]
      @p.run_steps(input).size.should == 1
      @p.run_steps(input)[0].should == :step_1
    end

    it "should take specified steps as input" do
      input = {:input_1 => "in", :steps => "step_2"}
      @p.run_steps(input).size.should == 1
      @p.run_steps(input)[0].should == :step_2
      input = {:input_1 => "in", :steps => ["step_1","step_2"]}
      @p.run_steps(input).size.should == 2
    end
  end
  describe "simple pipeline" do
    before(:each) do
      @simple = SimplePipeline.pipeline
    end
    it "should have one step" do
      @simple.steps.size.should == 1
    end
    it "should have nil as default steps" do
      @simple.default_steps.should == nil
    end
    it "should be prepared to run all steps" do
      @simple.run_steps.size.should == 1
      @simple.run_steps[0].should == :step_1
      inputs = {:steps => []}
      @simple.run_steps(inputs).size.should == 0
      inputs = {:steps => [:step_1]}
      @simple.run_steps.size.should == 1
      @simple.run_steps[0].should == :step_1
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
    it "should not run its step if no steps are to be run" do
      inputs = {:steps => [], :input_1 => 'input'}
      result = @simple.run inputs
      result.should == []
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
    it "should have correct output when all steps are run" do
      inputs = {:input_1 => 'aaa'}
      result = @two_step.run inputs
      result.should == ["step_1_aaa","step_2_aaa_out"]
    end
    it "should still have output from skipped steps" do
      inputs = {:steps => "step_2", :input_1 => 'aaa'}
      result = @two_step.run inputs
      result.should == ["step_2_aaa_out"]
    end
  end

  describe "pipeline with options" do
    before(:each) do
      @options_pipeline = OptionsPipeline.pipeline
    end
    it "should get options from options section" do
      parser = @options_pipeline.options_parser
      args = ["--input", "da/da/da"]
      @options_pipeline.parse_input args
      @options_pipeline.options[:input].should == "da/da/da"
    end
  end
end
