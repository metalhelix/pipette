

class String
  def append_filename appender
    self.split(".")[0..-2].join(".") + appender
  end
end
