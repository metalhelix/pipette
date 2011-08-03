
class SimplePipeline < Pipeline
  step :step_1 do
    input :input_1
    output :output_1 do |input| "#{input[:input_1]}_out" end
    run do |inputs, outputs|
      "run_#{inputs[:input_1]}_#{outputs[:output_1]}"
    end
  end
end

class TwoStepPipeline < Pipeline
  default_steps :step_1, :step_2
  step :step_1 do
    input :input_1
    output :output_1 do |input| "#{input[:input_1]}_out" end
    run do |inputs, outputs|
      "step_1_#{inputs[:input_1]}"
    end
  end
  step :step_2 do
    input :output_1
    output :output_2 do |input| "#{input[:output_1]}_out" end
    run do |inputs, outputs|
      "step_2_#{inputs[:output_1]}"
    end
  end
end
