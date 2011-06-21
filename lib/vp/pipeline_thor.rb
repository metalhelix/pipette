
require 'rubygems'
require 'thor'

require 'vp/gatk'
require 'vp/annotator'
require 'vp/vcf_filter'

module Vp
  class Pipeline < Thor

    class_option :yaml, :type => :string, :alias => "-y"

  end
end
