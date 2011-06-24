require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe Pipeline do
  describe "creation" do
    it "should initialize" do
      p = Pipeline.new
      p.steps.should == []
    end
  end
  describe "with steps" do
    before(:each) do 
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
  end
end
