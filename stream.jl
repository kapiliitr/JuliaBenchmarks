##!/nethome/kagarwal39/julia-0.3.1/julia/julia

addprocs(4);
NTIMES=3;
STREAM_ARRAY_SIZE = 5000000

isdefined(:STREAM_ARRAY_SIZE) || (STREAM_ARRAY_SIZE =	10000000);

if isdefined(:NTIMES)
    if NTIMES<=1
        NTIMES	= 10;
    end
end

isdefined(:NTIMES) || (NTIMES =	10);

isdefined(:OFFSET) || (OFFSET =	0);

HLINE = "-------------------------------------------------------------";

isdefined(:MIN) || (MIN(x,y)=((x)<(y)?(x):(y)));
isdefined(:MAX) || (MAX(x,y)=((x)>(y)?(x):(y)));

isdefined(:STREAM_TYPE) || (STREAM_TYPE = Float64);

a = Array(STREAM_TYPE,STREAM_ARRAY_SIZE+OFFSET);
b = Array(STREAM_TYPE,STREAM_ARRAY_SIZE+OFFSET);
c = Array(STREAM_TYPE,STREAM_ARRAY_SIZE+OFFSET);

avgtime = fill(0.0,4);
maxtime = fill(0.0,4);
mintime = {realmax(Float64),realmax(Float64),realmax(Float64),realmax(Float64)};

label = {"Copy:      ", "Scale:     ","Add:       ", "Triad:     "};

bytes = {
    2 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
    2 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
    3 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
    3 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE
    };

function main()
#    --- SETUP --- determine precision and check timing ---

    println(HLINE);
    println("Julia port of STREAM");
    println(HLINE);
    BytesPerWord = sizeof(STREAM_TYPE);
    println("This system uses ",BytesPerWord," bytes per array element.");

    println(HLINE);
    
    if isdefined(:N)
        println("*****  WARNING: ******");
        println("      It appears that you set the preprocessor variable N when compiling this code.");
        println("      This version of the code uses the preprocesor variable STREAM_ARRAY_SIZE to control the array size");
        println("      Reverting to default value of STREAM_ARRAY_SIZE=",STREAM_ARRAY_SIZE);
        println("*****  WARNING: ******");
    end

    println("Array size = ",STREAM_ARRAY_SIZE," (elements), Offset = ",OFFSET," (elements)");
    println("Memory per array = ",BytesPerWord * (STREAM_ARRAY_SIZE / 1024.0/1024.0)," MiB (= ",BytesPerWord * (STREAM_ARRAY_SIZE / 1024.0/1024.0/1024.0)," GiB)");
    println("Total memory required = ",(3.0 * BytesPerWord) * (STREAM_ARRAY_SIZE / 1024.0/1024.)," MiB (= ",(3.0 * BytesPerWord) * (STREAM_ARRAY_SIZE / 1024.0/1024.0/1024.0)," GiB)");
    println("Each kernel will be executed ",NTIMES," times.");
    println(" The *best* time for each kernel (excluding the first iteration)"); 
    println(" will be used to compute the reported bandwidth.");

#   Get initial value for system clock.
    @parallel for j=1:STREAM_ARRAY_SIZE
        a[j] = 1.0;
        b[j] = 2.0;
        c[j] = 0.0;
    end

    println(HLINE);

    if (quantum = checktick()) >= 1 
      	println("Your clock granularity/precision appears to be ",quantum," microseconds.");
    else
      	println("Your clock granularity appears to be less than one microsecond.");
      	quantum = 1;
    end

    t = @elapsed @parallel for j = 1:STREAM_ARRAY_SIZE
    		a[j] = 2.0E0 * a[j];
    end
    t = 1.0E6 * t;    

    println("Each test below will take on the order of ",convert(Int64,round(t))," microseconds.");
    println("   (= ",convert(Int64,round(t/quantum))," clock ticks)");
    println("Increase the size of the arrays if this shows that");
    println("you are not getting at least 20 clock ticks per test.");

    println(HLINE);

    println("WARNING -- The above is only a rough guideline.");
    println("For best results, please be sure you know the");
    println("precision of your system timer.");
    println(HLINE);
    
#   --- MAIN LOOP --- repeat test cases NTIMES times ---

    times = Array(Float64,(4,NTIMES))
    scalar = 3.0;
 
    for k=1:NTIMES
    
        if isdefined(:TUNED)
            times[1,k] = @elapsed tuned_STREAM_Copy();
        else
            times[1,k] = @elapsed @parallel for j=1:STREAM_ARRAY_SIZE
                              c[j] = a[j];
                          end
        end
      
        if isdefined(:TUNED)
            times[2,k] = @elapsed tuned_STREAM_Scale(scalar);
        else
            times[2,k] = @elapsed @parallel for j=1:STREAM_ARRAY_SIZE
                              b[j] = scalar*c[j];
                          end
        end        

        if isdefined(:TUNED)
            times[3,k] = @elapsed tuned_STREAM_Add();
        else
            times[3,k] = @elapsed @parallel for j=1:STREAM_ARRAY_SIZE
                              c[j] = a[j]+b[j];
                          end
        end

        if isdefined(:TUNED)
            times[4,k] = @elapsed tuned_STREAM_Triad(scalar);
        else
            times[4,k] = @elapsed @parallel for j=1:STREAM_ARRAY_SIZE
                               a[j] = b[j]+scalar*c[j];
                          end
        end
    end	

#   --- SUMMARY ---

    for k=2:NTIMES # note -- skip first iteration
	      for j=1:4
	          avgtime[j] = avgtime[j] + times[j,k];
      	    mintime[j] = MIN(mintime[j], times[j,k]);
      	    maxtime[j] = MAX(maxtime[j], times[j,k]);
        end
    end

    println("Function    Best Rate MB/s  Avg time     Min time     Max time");
    for j=1:4
    		avgtime[j] = avgtime[j]/(NTIMES-1);
#  	  	println(1.0E-06 * bytes[j]/mintime[j],"  ",avgtime[j],"  ",mintime[j],"  ",maxtime[j],"\n",label[j]);
        @printf("%s%12.1f  %11.6f  %11.6f  %11.6f\n", label[j],1.0E-06 * bytes[j]/mintime[j],avgtime[j],mintime[j],maxtime[j]);
    end

    println(HLINE);

#   --- Check Results ---
    checkSTREAMresults();
    println(HLINE);

    return 0;
end

M	= 20

function checktick()
    timesfound = Array(Float64,M);
    t2 = realmax(Float64);

#   Collect a sequence of M unique time values from the system.

    for i = 1:M
      	while (t2=@elapsed ()) < 1.0E-6
        end
      	timesfound[i] = t2;
    end

#   Determine the minimum difference between these M values.
#   This result will be our estimate (in microseconds) for the
#   clock granularity.
    
    minDelta = timesfound[1];
    for i = 2:M
      	minDelta = MIN(minDelta, timesfound[i]);
    end

   return convert(Int64, round(minDelta*1.0E6));
end

function checkSTREAMresults()
    
# reproduce initialization
  	aj = 1.0;
  	bj = 2.0;
  	cj = 0.0;
# a[] is modified during timing check
  	aj = 2.0E0 * aj;
# now execute timing loop
  	scalar = 3.0;
	  for k=1:NTIMES
        cj = aj;
        bj = scalar*cj;
        cj = aj+bj;
        aj = bj+scalar*cj;
    end

# accumulate deltas between observed and expected results
  	aSumErr = 0.0;
  	bSumErr = 0.0;
  	cSumErr = 0.0;
	  for j=1:STREAM_ARRAY_SIZE
    		aSumErr += abs(a[j] - aj);
		    bSumErr += abs(b[j] - bj);
    		cSumErr += abs(c[j] - cj);
	  end

    aAvgErr = aSumErr / STREAM_ARRAY_SIZE;
  	bAvgErr = bSumErr / STREAM_ARRAY_SIZE;
  	cAvgErr = cSumErr / STREAM_ARRAY_SIZE;

    epsilon = 1.e-6;
  	if sizeof(STREAM_TYPE) == 4
    		epsilon = 1.e-6;
		elseif sizeof(STREAM_TYPE) == 8
  	  	epsilon = 1.e-13;
	  else
  		  println("WEIRD: sizeof(STREAM_TYPE) = %lu\n",sizeof(STREAM_TYPE));
    end

  	err = 0;
  	if abs(aAvgErr/aj) > epsilon
    		err++;
		    println ("Failed Validation on array a[], AvgRelAbsErr > epsilon (",epsilon,")");
    		println ("     Expected Value: ",aj,", AvgAbsErr: ",aAvgErr,", AvgRelAbsErr: ",abs(aAvgErr)/aj);
    		ierr = 0;
    		for j=1:STREAM_ARRAY_SIZE
      			if abs(a[j]/aj-1.0) > epsilon
			          ierr++;
                if isdefined(:VERBOSE)
            				if ierr < 10
	  	  		            println("         array a: index: ",j,", expected: ",aj,", observed: ",a[j],", relative error: ",abs((aj-a[j])/aAvgErr));
            				end
                end
      			end
    		end
    		println("     For array a[], ",ierr," errors were found.");
    end

    if abs(bAvgErr/bj) > epsilon
    		err++;
		    println ("Failed Validation on array b[], AvgRelAbsErr > epsilon (",epsilon,")");
    		println ("     Expected Value: ",bj,", AvgAbsErr: ",bAvgErr,", AvgRelAbsErr: ",abs(bAvgErr)/bj);
    		println ("     AvgRelAbsErr > Epsilon (",epsilon,")");
    		ierr = 0;
		    for j=1:STREAM_ARRAY_SIZE
			      if abs(b[j]/bj-1.0) > epsilon
        				ierr++;
                if isdefined(:VERBOSE)
				            if ierr < 10
              					println("         array b: index: ",j,", expected: ",bj,", observed: ",b[j],", relative error: ",abs((bj-b[j])/bAvgErr));
				            end
                end
			      end
		    end
    		println("     For array b[], ",ierr," errors were found.");
	  end

	  if abs(cAvgErr/cj) > epsilon
    		err++;
		    println ("Failed Validation on array c[], AvgRelAbsErr > epsilon (%e)\n",epsilon);
    		println ("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n",cj,cAvgErr,abs(cAvgErr)/cj);
    		println ("     AvgRelAbsErr > Epsilon (%e)\n",epsilon);
    		ierr = 0;
    		for j=1:STREAM_ARRAY_SIZE
      			if abs(c[j]/cj-1.0) > epsilon
        				ierr++;
                if isdefined(:VERBOSE)
            				if ierr < 10
              					println("         array c: index: ",j,", expected: ",cj,", observed: ",c[j],", relative error: ",abs((cj-c[j])/cAvgErr));
				            end
                end
			      end
		    end
    		println("     For array c[], ",ierr," errors were found.");
    end

  	if err == 0
	    	println ("Solution Validates: avg error less than ",epsilon," on all three arrays");
  	end

    if isdefined(:VERBOSE)
      	println ("Results Validation Verbose Results: ");
      	println ("    Expected a(1), b(1), c(1): ",aj," ",bj," ",cj," ");
	      println ("    Observed a(1), b(1), c(1): ",a[i]," ",b[1]," ",c[1]," ");
      	println ("    Rel Errors on a, b, c:     ",abs(aAvgErr/aj)," ",abs(bAvgErr/bj)," ",abs(cAvgErr/cj)," ");
    end
end

# stubs for "tuned" versions of the kernels

isdefined(:TUNED) ||
function tuned_STREAM_Copy()
    @parallel for j=1:STREAM_ARRAY_SIZE
        c[j] = a[j];
    end
end

isdefined(:TUNED) ||
function tuned_STREAM_Scale(scalar::STREAM_TYPE)
	  @parallel for j=1:STREAM_ARRAY_SIZE
  	    b[j] = scalar*c[j];
    end
end

isdefined(:TUNED) ||
function tuned_STREAM_Add()
	  @parallel for j=1:STREAM_ARRAY_SIZE
  	    c[j] = a[j]+b[j];
    end
end

isdefined(:TUNED) ||
function tuned_STREAM_Triad(scalar::STREAM_TYPE)
	  @parallel for j=1:STREAM_ARRAY_SIZE
  	    a[j] = b[j]+scalar*c[j];
    end
end

# end of stubs for the "tuned" versions of the kernels

# Calling the main function
main()
