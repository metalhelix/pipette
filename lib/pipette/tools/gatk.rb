class Hash
  def to_options
    self.to_a.inject("") do |param_string, params| 
      if params[1].kind_of? Array
        params[1].each {|same_opt| param_string << " " << params[0] << " " << same_opt}
      else
        param_string << " " << params.join(" ") 
      end
      param_string.strip
    end
  end
end

class Array
  def to_options
    self.join(" ")
  end
end

class GATK
  def initialize options = {}
    defaults = { :verbose => true,
                 :cores => 4,
                 :java_params => ["-Xmx5g"],
                 :gatk_params => {"-et" => "NO_ET"}
               }
    @options = options.merge defaults
    check_options @options
  end

  def can_handle_multithreading? tool
    single_thread_only = ["RealignerTargetCreator", "IndelRealigner", "TableRecalibration", "VariantFiltration"]
    !single_thread_only.include? tool
  end

  def requires_reference_genome? tool
    true
  end

  def check_options options
    raise "ERROR GATK JAR not found at:#{options[:gatk]}." unless File.exists? options[:gatk]
  end

  def execute gatk_parameters
    if !gatk_parameters["-T"]
      raise "No GATK Tool provided. Add -T option"
    end
    java_options = @options[:java_params] ? @options[:java_params].to_options : ""
    gatk_options = @options[:gatk_params] ? @options[:gatk_params].to_options : ""

    gatk_options << " -nt #{@options[:cores]}" if can_handle_multithreading? gatk_parameters["-T"]
    gatk_options << " -R #{@options[:reference]}" if requires_reference_genome? gatk_parameters["-T"]

    gatk_params_string = gatk_parameters.to_options

    gatk_jar = @options[:gatk]

    gatk_call = "java #{java_options} -jar #{gatk_jar} #{gatk_options} #{gatk_params_string}"
    puts gatk_call if @options[:verbose]
    result = system(gatk_call)
  end
end
