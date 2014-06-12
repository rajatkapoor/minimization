# = minimization.rb -
# Minimization- Minimization algorithms on pure Ruby
# Copyright (C) 2010 Claudio Bustos
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
# Algorithms for unidimensional minimization
require 'text-table'
module Minimization
  VERSION = "0.2.1"
  @@gsl_present = false
  FailedIteration = Class.new(Exception)

  gspec = Gem::Specification.find_all_by_name("gsl")[0];
  if gspec !=nil
    @@gsl_present = true
    @@gsl_version = gspec.version
  end
  
  def self.has_gsl? #:nodoc:
    @@gsl_present
  end


# Base class for unidimensional minimizers

  class Unidimensional
    # Default value for error on f(x)
    EPSILON = 1e-6
    # Default number of maximum iterations
    MAX_ITERATIONS = 100
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
      raise "first argument should be lower than second" if lower>=upper
      @lower = lower
      @upper = upper
      @proc = proc
      golden = 0.3819660;
      @expected = @lower + golden * (@upper - @lower);
      @max_iteration = MAX_ITERATIONS
      @epsilon = EPSILON
      @iterations = 0
      @log = []
      @log_header = %w{I xl xh f(xl) f(xh) dx df(x)}
      
    end
    # Set expected value
    def expected=(v)
      @expected = v
    end
    def log_summary
      @log.join("\n")
    end
    # Convenience method to minimize
    # == Parameters:
    # * <tt>lower</tt>: Lower possible value
    # * <tt>upper</tt>: Higher possible value
    # * <tt>expected</tt>: Optional expected value. Faster the search is near correct value.
    # * <tt>&block</tt>: Block with function to minimize
    # == Usage:
    #   minimizer = Minimization::GoldenSection.minimize(-1000, 1000) {|x|
    #             x**2 }
    # 
    def self.minimize(lower, upper, expected = nil, &block)
      minimizer = new(lower, upper, block)
      minimizer.expected = expected unless expected.nil?
      raise FailedIteration unless minimizer.iterate
      minimizer
    end
    # Iterate to find the minimum
    def iterate
      raise "You should implement this"
    end
    def f(x)
      @proc.call(x)
    end
  end
  # Classic Newton-Raphson minimization method.  
  # Requires first and second derivative
  # == Usage
  #   f   = lambda {|x| x**2}
  #   fd  = lambda {|x| 2x}
  #   fdd = lambda {|x| 2}
  #   min = Minimization::NewtonRaphson.new(-1000, 1000, f, fd, fdd)
  #   min.iterate
  #   min.x_minimum
  #   min.f_minimum
  #   
  class NewtonRaphson < Unidimensional
    # == Parameters:
    # * <tt>lower</tt>: Lower possible value
    # * <tt>upper</tt>: Higher possible value
    # * <tt>proc</tt>: Original function
    # * <tt>proc_1d</tt>: First derivative
    # * <tt>proc_2d</tt>: Second derivative
    # 
    def initialize(lower, upper, proc, proc_1d, proc_2d)
      super(lower, upper, proc)
      @proc_1d = proc_1d
      @proc_2d = proc_2d
    end
    # Raises an error
    def self.minimize(*args)
      raise "You should use #new and #iterate"
    end
    def iterate
      # First
      x_prev = @lower
      x = @expected
      failed = true
      k = 0
      while (x-x_prev).abs > @epsilon and k<@max_iteration
        k+=1
        x_prev = x
        x = x-(@proc_1d.call(x).quo(@proc_2d.call(x)))
        f_prev = f(x_prev)
        f = f(x)
        x_min, x_max = [x, x_prev].min,  [x, x_prev].max 
        f_min, f_max = [f, f_prev].min,  [f, f_prev].max 
        @log << [k, x_min, x_max, f_min, f_max, (x_prev-x).abs, (f-f_prev).abs]
      end
      raise FailedIteration, "Not converged" if k>=@max_iteration
      @x_minimum = x;
      @f_minimum = f(x);
    end
  end
  # = Golden Section Minimizer.
  # Basic minimization algorithm. Slow, but robust.
  # See Unidimensional for methods.
  # == Usage.
  #  require 'minimization'
  #  min = Minimization::GoldenSection.new(-1000, 20000, proc {|x| (x+1)**2}
  #  min.expected=1.5  # Expected value
  #  min.iterate
  #  min.x_minimum
  #  min.f_minimum
  #  min.log
  class GoldenSection < Unidimensional
    # Start the iteration
    def gsl_iterate
    	@Gfun = GSL::Function.alloc(@proc)
		  k = 0
      @Gmin = FMinimizer.alloc(FMinimizer::GOLDENSECTION)
      m = @expected
      a = @lower
      b = @upper
      @Gmin.set(@Gfun,m,a,b)
      begin
      k += 1
      status = @Gmin.iterate
      status = @Gmin.test_interval(@epsilon, 0.0)
      puts("Converged:") if status == GSL::SUCCESS
      a = @Gmin.x_lower
      b = @Gmin.x_upper
      m = @Gmin.x_minimum
      end while status == GSL::CONTINUE and k < @max_iteration
    end

    def iterate
      ax = @lower
      bx = @expected
      cx = @upper
      c  =  (3-Math::sqrt(5)).quo(2);
      r  =  1-c;

      x0 = ax;
      x3 = cx;
      if ((cx-bx).abs > (bx-ax).abs)
        x1 = bx;
        x2 = bx + c*(cx-bx);
      else
        x2 = bx;
        x1 = bx - c*(bx-ax);
      end
      f1 = f(x1);
      f2 = f(x2);

      k = 1;
      while (x3-x0).abs > @epsilon and k<@max_iteration
        if f2 < f1
          x0 = x1;
          x1 = x2;
          x2 = r*x1 + c*x3;   # x2 = x1+c*(x3-x1)
          f1 = f2;
          f2 = f(x2);
        else
          x3 = x2;
          x2 = x1;
          x1 = r*x2 + c*x0;   # x1 = x2+c*(x0-x2)
          f2 = f1;
          f1 = f(x1);
        end
        @log << [k, x3, x0, f1, f2, (x3-x0).abs, (f1-f2).abs]
        
        k +=1;
      end

      if f1 < f2
        @x_minimum = x1;
        @f_minimum = f1;
      else
        @x_minimum = x2;
        @f_minimum = f2;
      end
      true
    end

  end

  # Direct port of Brent algorithm found on GSL.
  # See Unidimensional for methods.
  # == Usage
  #  min = Minimization::Brent.new(-1000, 20000, proc {|x| (x+1)**2}
  #  min.expected=1.5  # Expected value
  #  min.iterate
  #  min.x_minimum
  #  min.f_minimum
  #  min.log

  class Brent < Unidimensional
    GSL_SQRT_DBL_EPSILON=1.4901161193847656e-08
    def initialize(lower, upper,  proc)
      super

      @do_bracketing=true

      # Init

      golden = 0.3819660;      #golden = (3 - sqrt(5))/2

      v = @lower + golden * (@upper - @lower);
      w = v;

      @x_minimum = v ;
      @f_minimum = f(v) ;
      @x_lower = @lower
      @x_upper = @upper
      @f_lower = f(@lower) ;
      @f_upper = f(@lower) ;

      @v = v;
      @w = w;

      @d = 0;
      @e = 0;
      @f_v = f(v)
      @f_w = @f_v
    end

    def expected=(v)
      @x_minimum = v
      @f_minimum = f(v)
      @do_bracketing = false
    end
    
    def bracketing
      eval_max = 10
      f_left = @f_lower;
      f_right = @f_upper;
      x_left = @x_lower;
      x_right = @x_upper;
      golden = 0.3819660;      # golden = (3 - sqrt(5))/2 */
      nb_eval = 0

      if (f_right >= f_left)
        x_center = (x_right - x_left) * golden + x_left;
        nb_eval+= 1;
        f_center = f(x_center)
      else
        x_center = x_right ;
        f_center = f_right ;
        x_right = (x_center - x_left).quo(golden) + x_left;
        nb_eval+= 1;
        f_right = f(x_right);
      end


      begin
        @log << ["B#{nb_eval}", x_left, x_right, f_left, f_right, (x_left-x_right).abs, (f_left-f_right).abs]
        if (f_center < f_left )
          if (f_center < f_right)
            @x_lower  =  x_left;
            @x_upper = x_right;
            @x_minimum = x_center;
            @f_lower = f_left;
            @f_upper = f_right;
            @f_minimum = f_center;
            return true;
          elsif (f_center > f_right)
            x_left = x_center;
            f_left = f_center;
            x_center = x_right;
            f_center = f_right;
            x_right = (x_center - x_left).quo(golden) + x_left;
            nb_eval+=1;
            f_right = f(x_right);
          else # f_center == f_right */
            x_right = x_center;
            f_right = f_center;
            x_center = (x_right - x_left).quo(golden) + x_left;
            nb_eval+=1;
            f_center = f(x_center);
          end
        else # f_center >= f_left */
          x_right = x_center;
          f_right = f_center;
          x_center = (x_right - x_left) * golden + x_left;
          nb_eval+=1;
          f_center = f(x_center);
        end
      end while ((nb_eval < eval_max) and
      ((x_right - x_left) > GSL_SQRT_DBL_EPSILON * ( (x_right + x_left) * 0.5 ) + GSL_SQRT_DBL_EPSILON))
      @x_lower = x_left;
      @x_upper = x_right;
      @x_minimum = x_center;
      @f_lower = f_left;
      @f_upper = f_right;
      @f_minimum = f_center;
      return false;

    end
    # Start the minimization process
    # If you want to control manually the process, use brent_iterate
    def iterate
      k = 0
      bracketing if @do_bracketing
      while k<@max_iteration and (@x_lower-@x_upper).abs>@epsilon
        k+=1
        result = brent_iterate
        raise FailedIteration, "Error on iteration" if !result
        begin 
          @log << [k, @x_lower, @x_upper, @f_lower, @f_upper, (@x_lower-@x_upper).abs, (@f_lower-@f_upper).abs]
        rescue =>@e
          @log << [k, @e.to_s, nil, nil, nil, nil, nil]
        end
      end
      @iterations = k
      return true
    end
    # Generate one iteration.
    def brent_iterate
      x_left = @x_lower;
      x_right = @x_upper;

      z = @x_minimum;
      d = @e;
      e = @d;
      v = @v;
      w = @w;
      f_v = @f_v;
      f_w = @f_w;
      f_z = @f_minimum;

      golden = 0.3819660;      # golden = (3 - sqrt(5))/2 */

      w_lower = (z - x_left)
      w_upper = (x_right - z)

      tolerance =  GSL_SQRT_DBL_EPSILON * z.abs

      midpoint = 0.5 * (x_left + x_right)
      _p, q, r = 0, 0, 0
      if (e.abs > tolerance)

        # fit parabola */

        r = (z - w) * (f_z - f_v);
        q = (z - v) * (f_z - f_w);
        _p = (z - v) * q - (z - w) * r;
        q = 2 * (q - r);

        if (q > 0)
          _p = -_p
        else
          q = -q;
        end
        r = e;
        e = d;
      end

      if (_p.abs < (0.5 * q * r).abs and _p < q * w_lower and _p < q * w_upper)
        t2 = 2 * tolerance ;

        d = _p.quo(q);
        u = z + d;

        if ((u - x_left) < t2 or (x_right - u) < t2)
          d = (z < midpoint) ? tolerance : -tolerance ;
        end
      else

        e = (z < midpoint) ? x_right - z : -(z - x_left) ;
        d = golden * e;
      end

      if ( d.abs >= tolerance)
        u = z + d;
      else
        u = z + ((d > 0) ? tolerance : -tolerance) ;
      end

      @e = e;
      @d = d;

      f_u = f(u)

      if (f_u <= f_z)
        if (u < z)
          @x_upper = z;
          @f_upper = f_z;
        else
          @x_lower = z;
          @f_lower = f_z;
        end
        @v = w;
        @f_v = f_w;
        @w = z;
        @f_w = f_z;
        @x_minimum = u;
        @f_minimum = f_u;
        return true;
      else
        if (u < z)
          @x_lower = u;
          @f_lower = f_u;
          return true;
        else
          @x_upper = u;
          @f_upper = f_u;
          return true;
        end

        if (f_u <= f_w or w == z)
          @v = w;
          @f_v = f_w;
          @w = u;
          @f_w = f_u;
          return true;
        elsif f_u <= f_v or v == z or v == w
          @v = u;
          @f_v = f_u;
          return true;
        end
      end
      return false
    end

    def gsl_iterate
      @Gfun = GSL::Function.alloc(@proc)
      k = 0
      @Gmin = FMinimizer.alloc(FMinimizer::BRENT)
      m = @expected
      a = @lower
      b = @upper
      @Gmin.set(@Gfun,m,a,b)
      begin
      k += 1
      status = @Gmin.iterate
      status = @Gmin.test_interval(@epsilon, 0.0)
      puts("Converged:") if status == GSL::SUCCESS
      a = @Gmin.x_lower
      b = @Gmin.x_upper
      m = @Gmin.x_minimum
      end while status == GSL::CONTINUE and k < @max_iteration
    end

  end

  # = Bisection Method for Minimization.
  # See Unidimensional for methods.
  # == Usage.
  #  require 'minimization'
  #  min=Minimization::Bisection.new(1, 2  , proc {|x| (x)**3-(x)-2}
  #  min.iterate
  #  min.x_minimum
  #  min.f_minimum
  #  min.log
  # Source:
  #   * R.L. Burden, J. Faires: Numerical Analysis 
  class Bisection < Unidimensional

    def iterate()
      ax = @lower
      cx = @upper
      k = 0;
      while (ax-cx).abs > @epsilon and k<@max_iteration
        bx = (ax + cx).quo(2);
        fa = f(ax);
        fb = f(bx);
        fc = f(cx);
        if (fa*fb <0)
          cx = bx;
        else (fb*fc <0)
          ax = bx;
        end
        k +=1;
        @log << [k, ax.to_f, cx.to_f, f(ax).to_f, f(cx).to_f, (ax-cx).abs.to_f, (f(ax)-f(cx)).abs.to_f]
      end
      
      if (fa<fc)
        @x_minimum, @f_minimum = ax.to_f, f(ax).to_f;
      else 
        @x_minimum, @f_minimum = cx.to_f, f(cx).to_f;
      end
    end
  end

  #Quad Golden method
  class QuadGolden < Unidimensional
      REL_ERR_VAL = 1.0e-06;
      GSL_DBL_EPSILON = 2.2204460492503131e-16;
      ABS_ERR_VAL = 1.0e-10;
      GOLDEN_MEAN = 0.3819660112501052;
      GOLDEN_RATIO = 1.6180339887498950;
      golden = 0.3819660;
  
#    def intialize(lower, upper, proc)
#      super
# #     @x_prev_small = @expected;
#      @x_small = @expected;
#
#      @f_prev_small = f(@expected);
#      @f_small = f(@expected);
#
 #     @step_size = 0.0;
#      @stored_step = 0.0;
 #     @prev_stored_step = 0.0;
#      @num_iter = 0;
#
#      @x_lower = lower ;
#      @x_upper = upper ;
#      @f_lower = f(@x_lower) ;
#      @f_upper = f(@x_upper) ;
#      f = 0;
#      puts (lower, upper)
#    end

    def iterate()

      @x_prev_small = @expected;
      @x_small = @expected;

      @f_prev_small = f(@expected);
      @f_small = f(@expected);

      @step_size = 0.0;
      @stored_step = 0.0;
      @prev_stored_step = 0.0;
      @num_iter = 0;

      @x_lower = @lower ;
      @x_upper = @upper ;
      @f_lower = f(@x_lower) ;
      @f_upper = f(@x_upper) ;
      f = 0;

      x_m = @expected;
      f_m = f(@expected);

      x_l = @lower;
      x_u = @upper;

      x_small = @x_small;
      f_small = @f_small;

      x_prev_small = @x_prev_small;
      f_prev_small = @f_prev_small;
  
      stored_step = @stored_step; # update on exit 
      prev_stored_step = @prev_stored_step; #update on exit 
      step_size = @step_size; # update on exit 

      quad_step_size = prev_stored_step;
   
      x_midpoint = 0.5 * (x_l + x_u);
      tol = REL_ERR_VAL * x_m.abs + ABS_ERR_VAL; # total error tolerance 



      if (stored_step.abs - tol > -2.0 * GSL_DBL_EPSILON)
        #Fit quadratic 
        double c3 = (x_m - x_small) * (f_m - f_prev_small);
        double c2 = (x_m - x_prev_small) * (f_m - f_small);
        double c1 = (x_m - x_prev_small) * c2 - (x_m - x_small) * c3;

        c2 = 2.0 * (c2 - c3);

        if (c2.abs > GSL_DBL_EPSILON) # if( c2 != 0 ) 
          if (c2 > 0.0)
            c1 = -c1;
          end

          c2 =  c2.abs;

          quad_step_size = c1.quo(c2);
        else
    
        # Handle case where c2 ~=~ 0  */
        # Insure that the line search will NOT take a quadratic interpolation step in this iteration
        quad_step_size = stored_step;
        end

        prev_stored_step = stored_step;
        stored_step = step_size;
      end

      x_trial = x_m + quad_step_size;

      if (quad_step_size.abs <  (0.5 * prev_stored_step).abs && x_trial > x_l && x_trial < x_u)
      
        #/* Take quadratic interpolation step */
        step_size = quad_step_size;

        #/* Do not evaluate function too close to x_l or x_u */
        if ((x_trial - x_l) < 2.0 * tol || (x_u - x_trial) < 2.0 * tol)
          step_size = (x_midpoint >= x_m ? +1.0 : -1.0) * tol.abs;
        end

        #DEBUG_PRINTF(("quadratic step: %g\n", step_size));
        puts "quadratic step: " + step_size.to_s,"";
      
      elsif ((x_small != x_prev_small && x_small < x_m && x_prev_small < x_m) ||
             (x_small != x_prev_small && x_small > x_m && x_prev_small > x_m))
        #/* Take safeguarded function comparison step */
        #double outside_interval, inside_interval;
        if (x_small < x_m)
    
          outside_interval = x_l - x_m;
          inside_interval = x_u - x_m;
        else
    
          outside_interval = x_u - x_m;
          inside_interval = x_l - x_m;
        end
        if (inside_interval.abs <= tol)
            #/* Swap inside and outside intervals */
            double tmp = outside_interval;
            outside_interval = inside_interval;
            inside_interval = tmp;
        end

        #{
          double step = inside_interval;
          double scale_factor;

          if (outside_interval.abs < inside_interval.abs)
              scale_factor = 0.5 * Math::sqrt(-outside_interval.quo(inside_interval));
          else
              scale_factor = (5.0 / 11.0) * (0.1 - inside_interval.quo(outside_interval));
          end

          @stored_step = step;
          step_size = scale_factor * step;
        #}

        #DEBUG_PRINTF(("safeguard step: %g\n", step_size));
        puts "safeguarded step: " + step_size,"";
    
      else
        #/* Take golden section step */
        #double step;

        if (x_m < x_midpoint)
            step = x_u - x_m;
        else
            step = x_l - x_m;
        end
        @stored_step = step;
        step_size = GOLDEN_MEAN * step;

        #DEBUG_PRINTF(("golden step: %g\n", step_size));
      end

      #/* Do not evaluate function too close to x_minimum */
      if (step_size.abs > tol)
        x_eval = x_m + step_size;
      else
        x_eval = x_m + (step_size >= 0 ? +1.0 : -1.0) * tol.abs;
      end
      #/* Evaluate function at the new point x_eval */
      f_eval = f(x_eval);

      #/* Update {x,f}_lower, {x,f}_upper, {x,f}_prev_small, {x,f}_small, and {x,f}_minimum */
      if (f_eval <= f_m)
          if (x_eval < x_m)
              @x_upper = x_m;
              @f_upper = f_m;     
          else
              @x_lower = x_m;
              @f_upper = f_m;
          end

          @x_prev_small = x_small;
          @f_prev_small = f_small;

          @x_small = x_m;
          @f_small = f_m;

          @x_minimum = x_eval;
          @f_minimum = f_eval;
      else
          if (x_eval < x_m)
              @x_lower = x_eval;
              @f_lower = f_eval;
          else
                
              @x_upper = x_eval;
              @f_upper = f_eval;
          end
          if (f_eval <= f_small ||  (x_small - x_m).abs < 2.0 * GSL_DBL_EPSILON)

              @x_prev_small = x_small;
              @f_prev_small = f_small;

              @x_small = x_eval;
              @f_small = f_eval;

          elsif (f_eval <= f_prev_small ||
             (x_prev_small - x_m).abs < 2.0 * GSL_DBL_EPSILON ||
             (x_prev_small - x_small).abs < 2.0 * GSL_DBL_EPSILON)
      
              @x_prev_small = x_eval;
              @f_prev_small = f_eval;
          end
      
      end
    


    #/* Update stored values for next iteration */
    @prev_stored_step = prev_stored_step;
    @step_size = step_size;
    @num_iter+=1; 

    puts @num_iter.to_s + "Final State" + x_l.to_s , x_m.to_s , x_u.to_s;
    end

    def gsl_iterate
      @Gfun = GSL::Function.alloc(@proc)
      k = 0
      @Gmin = FMinimizer.alloc(FMinimizer::QUADGOLDEN)
      m = @expected
      a = @lower
      b = @upper
      @Gmin.set(@Gfun,m,a,b)
      begin
      k += 1
      status = @Gmin.iterate
      status = @Gmin.test_interval(@epsilon, 0.0)
      puts("Converged:") if status == GSL::SUCCESS
      a = @Gmin.x_lower
      b = @Gmin.x_upper
      m = @Gmin.x_minimum
      end while status == GSL::CONTINUE and k < @max_iteration
    end
  end
end