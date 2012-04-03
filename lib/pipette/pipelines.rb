
# hack to get all the subclasses for a class
# found here:
# http://stackoverflow.com/questions/436159/how-to-get-all-subclasses
class Class
  def subclasses
    result = []
    ObjectSpace.each_object(Class) { |klass| result << klass if klass < self }
    result
  end
end

# if the $PIPETTE_DIR environmental variable is set,
# add it to the pipeline directories.
# This might need to be changed if we want to support
# multiple directories in $PIPETTE_DIR
EXTERNAL_PIPELINES_DIR = ENV['PIPETTE_DIR']

PIPELINES_DIRS = [EXTERNAL_PIPELINES_DIR, File.join(File.dirname(__FILE__), 'pipelines')].compact

PIPELINES_DIRS.each do |dir|
  Dir[File.join(dir, '*.rb')].each {|file| require file }
end

all_pipelines = {}
Pipeline.subclasses.each do |subclass|
  pipeline = subclass.pipeline
  all_pipelines[pipeline.name] = subclass
end

# set our ALL_PIPELINES constant to be used
# in pipette.rb
ALL_PIPELINES = all_pipelines

