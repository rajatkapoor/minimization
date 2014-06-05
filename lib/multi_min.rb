 require 'text-table'
 module Minimization
  class Multidimensional
    # Default value for error on function
    #EPSILON = 1e-6
    # Default number of maximum iterations
    #MAX_ITERATIONS = 100
    # Minimum value for x
    attr_reader :x_minimum
    # Minimum value for f(x)
    attr_reader :f_minimum
    # Log of iterations. Should be an array
    attr_reader :log
    # Name of fields of log
    attr_reader :log_header
    # Absolute error on x
    attr_accessor :epsilon
    # Expected value. Fast minimum finding if set
    attr_reader :expected
    # Numbers of iterations
    attr_reader :iterations
    # Create a new minimizer
    def initialize(lower, upper, proc)
      @lower = lower
      @upper = upper
      @proc = proc
      @nvars = @proc.arity
      @params = @proc.parameters
    end

    def f(*vars)
      @proc.call(vars)
    end
  end
end