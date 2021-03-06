
require 'pipette/tools/vcf'

class VCFFilter
  def self.filter input, output
    vcf = VCF.new input
    vcf.copy_if(output) do |vcf_hash| 
      rtn = vcf_hash["DP2"] > 10
      rtn &= vcf_hash["genotype"] != "isHet"
      # if rtn && vcf_hash["genotype"] == "isHet"
      #  val  = (vcf_hash["AD"][1].to_f / vcf_hash["AD"][0].to_f)
      #  rtn = val > 2.0
        #puts "low divide: #{val}" unless rtn
      # end
      rtn
    end
    output
  end
end
