#!/usr/bin/env ruby


FILTER_MATCH = File.expand_path(File.join(File.dirname(__FILE__), "..", "bin", "fast_filter_matching_vcfs.rb"))

FILTER = File.expand_path(File.join(File.dirname(__FILE__), "..", "bin", "filter_vcf.rb"))

ANNOTATE = File.expand_path(File.join(File.dirname(__FILE__), "..", "bin", "annotate_vcf.rb"))

ANNOTATE_CONFIG = File.expand_path(File.join(File.dirname(__FILE__), "variant_config.yml"))

LINKS = File.expand_path(File.join(File.dirname(__FILE__), "..", "bin", "add_igv_links.rb"))

INDEX_VCF = "java -Xmx1500m -Djava.awt.headless=true -jar /n/local/stage/igv_tools/current/igvtools.jar index"


starting_dir = ARGV[0]


control_name = "Control"



def filter_matching_vcfs vcf_files, control_name, starting_dir
  control_snp_vcf_filename = File.expand_path(File.join(starting_dir, control_name, "#{control_name}.snps.vcf"))
  control_indel_vcf_filename = File.expand_path(File.join(starting_dir, control_name, "#{control_name}.indels.vcf"))
  vcf_files.each do |vcf_file|
    puts vcf_file
    dir_name = File.dirname(vcf_file)
    Dir.chdir(dir_name) do
      vcf_basename = File.basename(vcf_file)

      filter_file = vcf_basename.split(".")[1] == "indels" ? control_indel_vcf_filename : control_snp_vcf_filename

      command = "#{FILTER_MATCH} #{vcf_basename} #{filter_file}"
      puts command
      system(command)
    end
  end
end

def filter_vcfs vcf_files
  vcf_files.each do |vcf_file|
    puts vcf_file
    dir_name = File.dirname(vcf_file)
    Dir.chdir(dir_name) do
      vcf_basename = File.basename(vcf_file)

      filter_name = File.basename(vcf_file, File.extname(vcf_file)) + ".filter.vcf"

      command = "#{FILTER} #{vcf_basename} #{filter_name}"
      puts command
      system(command)
    end
  end
end

def annotate_vcfs vcf_files
  vcf_files.each do |vcf_file|
    puts vcf_file
    dir_name = File.dirname(vcf_file)
    Dir.chdir(dir_name) do
      vcf_basename = File.basename(vcf_file)

      command = "#{ANNOTATE} -y #{ANNOTATE_CONFIG} -i #{vcf_basename}"
      puts command
      system(command)
    end
  end
end

def add_links tsv_files
  tsv_files.each do |tsv_file|
    puts tsv_file
    dir_name = File.dirname(tsv_file)
    Dir.chdir(dir_name) do
      basename = File.basename(tsv_file)
      command = "#{LINKS} #{basename}"

      puts command
      system(command)

    end
  end
end

def copy_results results_files, results_dir
  system("mkdir -p #{results_dir}")
  results_files.each do |result_file|
    command = "cp #{result_file} #{results_dir}"
    puts command
    system(command)
  end
end

def index_vcf_files vcf_files
  vcf_files.each do |vcf_file|
    Dir.chdir(File.dirname(vcf_file)) do
      command = "#{INDEX_VCF} #{File.basename(vcf_file)}"
      puts command
      system(command)
    end
  end
end

vcf_files = Dir.glob(File.join(starting_dir, "**", "*.vcf"))
filter_matching_vcfs vcf_files, control_name, starting_dir

vcf_files = Dir.glob(File.join(starting_dir, "**", "*.unique_from_#{control_name}.*.vcf"))
filter_vcfs(vcf_files)

vcf_files = Dir.glob(File.join(starting_dir, "**", "*.unique_from_#{control_name}.*.filter.vcf"))
annotate_vcfs vcf_files

output_dir = File.join("results", "vcf_files")
copy_results vcf_files, output_dir

vcf_files = Dir.glob(File.join(output_dir, "*.vcf"))
index_vcf_files vcf_files

tsv_files = Dir.glob(File.join(starting_dir, "**", "*filter*.collapse.txt"))
add_links(tsv_files)

result_files = Dir.glob(File.join(starting_dir, "**", "*filter*.links.xls"))
output_dir = File.join("results", "linked_files")
copy_results result_files, output_dir

output_dir = File.join("results", "summary_files")
result_files = Dir.glob(File.join(starting_dir, "**", "*filter*.snpeff.summary.html"))
copy_results result_files, output_dir


