include("iridis_launcher.jl")

if(length(ARGS)!=2)
	println("Usage : julia fft.jl <N> <P>")
	quit()
end

N=parseint(ARGS[1])
P=parseint(ARGS[2])

bind_iridis_procs(P)

if((ceil(log2(N))-log2(N))!=0 && (ceil(log2(nprocs()-1))-log2(nprocs()-1))!=0)
	println("Array size(N) and Number of Processes must be powers of two!")
	quit()
end

function getBit(x,i)
	z= convert(Uint64,x)>>convert(Uint64,i)
        y=z&1
        return y
end

function reverseBits(x,N)
	for i=0:(N/2)-1
        	if(getBit(x,i)!=getBit(x,N-1-i))
                	y=(convert(Uint64,1)<<convert(Uint64,i))|(convert(Uint64,1)<<convert(Uint64,(N-1-i)))
                   	x=convert(Uint64,x)$y
               	end
        end
        x
end

@everywhere function localComputation(D::DArray)
	A=fft(localpart(D))
	localpart(D)[1]=A[1]
        for i=length(localpart(D)):-1:2
        	localpart(D)[i]=A[length(localpart(D))-i+2]
        end

	
end

@everywhere function fetch_data(D::DArray)
        return localpart(D)
end


@everywhere nthroots(n::Integer) = [ cospi(2k/n)+sinpi(2k/n)im for k = 0:n-1 ]

@everywhere function parallel_fft(D::DArray, r, Aux::DArray)
	#Proc to fetch data from
	pos=(myid()-2)%(2^r)
	if(pos<2^(r-1))
		proc_num=myid()+2^(r-1)
	else
		proc_num=myid()-2^(r-1)
	end
	B=fetch(@spawnat proc_num fetch_data(D))
	
	#Creating the corresponding roots of unity
	start_index=((myid()-2)*length(localpart(D)))%convert(Integer,(length(localpart(D))*(2^r))) +1
	Roots=nthroots(convert(Integer,length(localpart(D))*(2^r)))[start_index:(start_index+length(localpart(D))-1)]
	
	if(proc_num>myid())
		for i = 1:length(localpart(D))
			localpart(Aux)[i]=localpart(D)[i]+Roots[i]*B[i]
		end	
	else
		for i = 1:length(localpart(D))
			localpart(Aux)[i]=B[i]+Roots[i]*localpart(D)[i]
		end
	end	
end

@everywhere function copy_array(D::DArray, Aux::DArray)
	for i =1:length(localpart(D))
		localpart(D)[i]=localpart(Aux)[i]
	end
end


function verify()
	seq_fft=fft(A)
	seq_fft=[seq_fft[1],reverse(seq_fft[2:end])]

	par_fft=convert(Array,D)

	if(length(par_fft)!=length(seq_fft))
		println("ERROR Detected : Size of arrays don't match")
	else
		flag=0
		for i = 1:length(par_fft)
			if(abs(par_fft[i]-seq_fft[i])>0.01)
				println("WRONG Result")
				println(seq_fft)
				println(par_fft)
				flag=1
				break
			end
		end

		if(flag==0)
			println("Result is correct!")
		end
	end

end

A=rand(N)
B=zeros(Complex,length(A))

for i=0:nprocs()-2
        X=Array((Number,Number),convert(Integer,N/(nprocs()-1)))
       	count=1
        for j=i*(N/(nprocs()-1)):(i+1)*(N/(nprocs()-1))-1
                X[count]=(reverseBits(j,log2(N)),convert(Complex,A[reverseBits(j,log2(N))+1]))
                count+=1
        end   
        Y=sort(X)
        count=1
        for j=i*(N/(nprocs()-1)):(i+1)*(N/(nprocs()-1))-1
              	B[j+1]=Y[count][2]
               	count+=1
        end
end

D=distribute(B)

#Local computation of chunks at each processor
time_1=@elapsed @sync {@spawnat p localComputation(D) for p in procs(D)}

Aux=zeros(Complex,N)
Aux=distribute(Aux)

# Communication between processors
time_2=@elapsed for r=1:convert(Integer,log2(nprocs()-1))

	@sync {@spawnat p parallel_fft(D,r,Aux) for p in procs(D)}

	@sync {@spawnat p copy_array(D,Aux) for p in procs(D)}

end

verify()

println("Time taken : ",time_1+time_2)
println("Number of processes : ",nprocs()-1)
