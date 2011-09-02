Dir[File.dirname(__FILE__) + '/../tools/*.rb'].each do |file|
  puts "requiring: #{file}"
  require File.basename(file, File.extname(file))
end
