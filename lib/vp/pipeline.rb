class Pipeline
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
