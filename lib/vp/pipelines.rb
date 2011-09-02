Dir[File.dirname(__FILE__) + '/../pipelines/*.rb'].each do |file|
  require File.basename(file, File.extname(file))
end
