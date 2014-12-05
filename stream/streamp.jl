#!/nethome/kagarwal39/julia-0.3.3/julia/julia

isdefined(:STREAM_ARRAY_SIZE) || (const STREAM_ARRAY_SIZE =	64000000);
isdefined(:NTIMES) || (const NTIMES =	10);
isdefined(:OFFSET) || (const OFFSET =	0);
isdefined(:STREAM_TYPE) || (const STREAM_TYPE = Float64);

#numargs = countnz(ARGS);
const PARALLEL = 1;
#addprocs(4);
NUMPROCS=parseint(ARGS[1]);
#println(NUMPROCS, " threads");
addprocs(NUMPROCS);


#if numargs>3
#    println("Invalid number of arguments. ./bfs [STREAM_ARRAY_SIZE NTIMES NUM_WORKERS]");
#    exit(-1);
#end
#if numargs>=1
#    const STREAM_ARRAY_SIZE = parseint(ARGS[1]);
#end
#if numargs>=2
#    const NTIMES = parseint(ARGS[2]);
#end
#if numargs==3
#    const PARALLEL = 1;
#    addprocs(parseint(ARGS[3]));
#end
#
#if isdefined(:NTIMES)
#    if NTIMES<=1
#        NTIMES	= 10;
#    end
#end

isdefined(:MIN) || (MIN(x,y)=((x)<(y)?(x):(y)));
isdefined(:MAX) || (MAX(x,y)=((x)>(y)?(x):(y)));

function main()

    HLINE = "-------------------------------------------------------------";

    a = fill(1.0,STREAM_ARRAY_SIZE::Int64+OFFSET::Int64);
    b = fill(2.0,STREAM_ARRAY_SIZE::Int64+OFFSET::Int64);
    c = fill(0.0,STREAM_ARRAY_SIZE::Int64+OFFSET::Int64);

    if isdefined(:PARALLEL)
        a = distribute(a);
        b = distribute(b);
        c = distribute(c);
    end

    avgtime = fill(0.0,4);
    maxtime = fill(0.0,4);
    mintime = {realmax(Float64),realmax(Float64),realmax(Float64),realmax(Float64)};

    label = {"Copy:      ", "Scale:     ","Add:       ", "Triad:     "};

    bytes = {
        2 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE::Int64,
        2 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE::Int64,
        3 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE::Int64,
        3 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE::Int64
        };

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
        println("      Reverting to default value of STREAM_ARRAY_SIZE=",STREAM_ARRAY_SIZE::Int64);
        println("*****  WARNING: ******");
    end

    println("Array size = ",STREAM_ARRAY_SIZE::Int64," (elements), Offset = ",OFFSET::Int64," (elements)");
    println("Memory per array = ",BytesPerWord * (STREAM_ARRAY_SIZE::Int64 / 1024.0/1024.0)," MiB (= ",BytesPerWord * (STREAM_ARRAY_SIZE::Int64 / 1024.0/1024.0/1024.0)," GiB)");
    println("Total memory required = ",(3.0 * BytesPerWord) * (STREAM_ARRAY_SIZE::Int64 / 1024.0/1024.)," MiB (= ",(3.0 * BytesPerWord) * (STREAM_ARRAY_SIZE::Int64 / 1024.0/1024.0/1024.0)," GiB)");
    println("Each kernel will be executed ",NTIMES::Int64," times.");
    println(" The *best* time for each kernel (excluding the first iteration)"); 
    println(" will be used to compute the reported bandwidth.");

    if mod(STREAM_ARRAY_SIZE::Int64,nworkers()) != 0
        println("*****  ERROR: ******");
        println("      The preprocessor variable STREAM_ARRAY_SIZE is not a multiple of the number of worker processes.");
        println("*****  ERROR: ******");
        exit(-1);
    end

    println(HLINE);

    if (quantum = checktick()) >= 1 
      	println("Your clock granularity/precision appears to be ",quantum," microseconds.");
    else
      	println("Your clock granularity appears to be less than one microsecond.");
      	quantum = 1;
    end

    if isdefined(:PARALLEL)
        t = @elapsed @sync { (@spawnat p double(localpart(a))) for p=procs(a) }
    else
        t = @elapsed for j = 1:STREAM_ARRAY_SIZE::Int64
        		a[j] = 2.0E0 * a[j];
        end
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

    times = Array(Float64,(4,NTIMES::Int64))
    scalar = 3.0;

    for k=1:NTIMES::Int64

        if isdefined(:PARALLEL)
            times[1,k] = @elapsed @sync { (@spawnat p parallel_STREAM_Copy(localpart(a),localpart(c))) for p=procs(a) }
        else
            times[1,k] = @elapsed for j=1:STREAM_ARRAY_SIZE::Int64
                              c[j] = a[j];
                          end
        end
      
        if isdefined(:PARALLEL)
            times[2,k] = @elapsed @sync { (@spawnat p parallel_STREAM_Scale(localpart(c),localpart(b),scalar)) for p=procs(c) }
        else
            times[2,k] = @elapsed for j=1:STREAM_ARRAY_SIZE::Int64
                              b[j] = scalar*c[j];
                          end
        end        

        if isdefined(:PARALLEL)
            times[3,k] = @elapsed @sync { (@spawnat p parallel_STREAM_Add(localpart(a),localpart(b),localpart(c))) for p=procs(a) }
        else
            times[3,k] = @elapsed for j=1:STREAM_ARRAY_SIZE::Int64
                              c[j] = a[j]+b[j];
                          end
        end

        if isdefined(:PARALLEL)
            times[4,k] = @elapsed @sync { (@spawnat p parallel_STREAM_Triad(localpart(b),localpart(c),localpart(a),scalar)) for p=procs(b) }
        else
            times[4,k] = @elapsed for j=1:STREAM_ARRAY_SIZE::Int64
                               a[j] = b[j]+scalar*c[j];
                          end
        end
    end	

#   --- SUMMARY ---

    for k=2:NTIMES::Int64 # note -- skip first iteration
	      for j=1:4
	          avgtime[j] = avgtime[j] + times[j,k];
      	    mintime[j] = MIN(mintime[j], times[j,k]);
      	    maxtime[j] = MAX(maxtime[j], times[j,k]);
        end
    end

    println("Function    Best Rate MB/s  Avg time     Min time     Max time");
    for j=1:4
    		avgtime[j] = avgtime[j]/(NTIMES::Int64-1);
        @printf("%s%12.1f  %11.6f  %11.6f  %11.6f\n", label[j],1.0E-06 * bytes[j]/mintime[j],avgtime[j],mintime[j],maxtime[j]);
    end

    println(HLINE);

#   --- Check Results ---
    checkSTREAMresults(a,b,c);
    println(HLINE);

#    exit(0);
end

function checktick()
    M	= 20
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

function checkSTREAMresults(a,b,c)

# reproduce initialization
  	aj = 1.0;
  	bj = 2.0;
  	cj = 0.0;
# a[] is modified during timing check
  	aj = 2.0E0 * aj;
# now execute timing loop
  	scalar = 3.0;
	  for k=1:NTIMES::Int64
        cj = aj;
        bj = scalar*cj;
        cj = aj+bj;
        aj = bj+scalar*cj;
    end

# accumulate deltas between observed and expected results
  	aSumErr = 0.0;
  	bSumErr = 0.0;
  	cSumErr = 0.0;
    
    if isdefined(:PARALLEL)
        aSumErr = reduce(+, map(fetch,{ (@spawnat p calcErr(localpart(a),aj)) for p=procs(a) }));
        bSumErr = reduce(+, map(fetch,{ (@spawnat p calcErr(localpart(b),bj)) for p=procs(b) }));
        cSumErr = reduce(+, map(fetch,{ (@spawnat p calcErr(localpart(c),cj)) for p=procs(c) }));
    else
    	  for j=1:STREAM_ARRAY_SIZE::Int64
        		aSumErr += abs(a[j] - aj);
		        bSumErr += abs(b[j] - bj);
        		cSumErr += abs(c[j] - cj);
    	  end
    end

    aAvgErr = aSumErr / STREAM_ARRAY_SIZE::Int64;
  	bAvgErr = bSumErr / STREAM_ARRAY_SIZE::Int64;
  	cAvgErr = cSumErr / STREAM_ARRAY_SIZE::Int64;

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
    		err += 1;
		    println ("Failed Validation on array a[], AvgRelAbsErr > epsilon (",epsilon,")");
    		println ("     Expected Value: ",aj,", AvgAbsErr: ",aAvgErr,", AvgRelAbsErr: ",abs(aAvgErr)/aj);
    		ierr = 0;
    		for j=1:STREAM_ARRAY_SIZE::Int64
      			if abs(a[j]/aj-1.0) > epsilon
			          ierr += 1;
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
    		err += 1;
		    println ("Failed Validation on array b[], AvgRelAbsErr > epsilon (",epsilon,")");
    		println ("     Expected Value: ",bj,", AvgAbsErr: ",bAvgErr,", AvgRelAbsErr: ",abs(bAvgErr)/bj);
    		println ("     AvgRelAbsErr > Epsilon (",epsilon,")");
    		ierr = 0;
		    for j=1:STREAM_ARRAY_SIZE::Int64
			      if abs(b[j]/bj-1.0) > epsilon
        				ierr += 1;
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
    		err += 1;
		    println ("Failed Validation on array c[], AvgRelAbsErr > epsilon (",epsilon,")");
    		println ("     Expected Value: ",cj,", AvgAbsErr: ",cAvgErr,", AvgRelAbsErr: ",abs(cAvgErr)/cj);
    		println ("     AvgRelAbsErr > Epsilon (",epsilon,")");
    		ierr = 0;
    		for j=1:STREAM_ARRAY_SIZE::Int64
      			if abs(c[j]/cj-1.0) > epsilon
        				ierr += 1;
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

# stubs for "parallel" versions of the kernels

~isdefined(:PARALLEL) ||
@everywhere function double(a::Array)
    for i=1:size(a)[1]
        a[i] = 2*a[i]
    end
end

~isdefined(:PARALLEL) ||
@everywhere function parallel_STREAM_Copy(a::Array,c::Array)
    for i=1:size(c)[1]
        c[i] = a[i]
    end
end

~isdefined(:PARALLEL) ||
@everywhere function parallel_STREAM_Scale(c::Array,b::Array,s::Float64)
    for i=1:size(b)[1]
  	    b[i] = s*c[i];
    end
end

~isdefined(:PARALLEL) ||
@everywhere function parallel_STREAM_Add(a::Array,b::Array,c::Array)
    for i=1:size(c)[1]
  	    c[i] = a[i]+b[i];
    end
end

~isdefined(:PARALLEL) ||
@everywhere function parallel_STREAM_Triad(b::Array,c::Array,a::Array,s::Float64)
    for i=1:size(a)[1]
  	    a[i] = b[i]+s*c[i];
    end
end

~isdefined(:PARALLEL) ||
@everywhere function calcErr(a::Array,aj::Float64)
    sum = 0.0
    for i=1:size(a)[1]
  	    sum += abs(a[i] - aj);
    end
    return sum;
end

# end of stubs for the "parallel" versions of the kernels

# Calling the main function
main()
