
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
  step :step_1 do
    run do |inputs, outputs|
    end
  end
  step :step_2 do
    run do |inputs, outputs|
    end
  end
end
