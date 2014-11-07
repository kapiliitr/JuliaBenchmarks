addprocs(4)

#Define variables on all procs
@everywhere LCG_MUL64=6364136223846793005
@everywhere LCG_ADD64=1
@everywhere TableSize=1024*4096
@everywhere logTableSize=log2(TableSize)
@everywhere numProcs=length(workers())
@everywhere logNumProcs=log2(numProcs)
@everywhere maxPendingUpdates=1024
@everywhere Buckets=Array(Uint64,numProcs,maxPendingUpdates)


#Declare a darray on master node
HPCC_Table=Array(Uint64,TableSize)


#Initialize the darray such that HPCC_Table[i]=i
for i in 1:TableSize
	HPCC_Table[i]=i-1
end
HPCC_Table=distribute(HPCC_Table)

# Utility routine to start random number generator at Nth step
@everywhere function HPCC_starts(n)
	global LCG_MUL64+=0
        global LCG_ADD64+=0
        mul_k=convert(Uint64,LCG_MUL64)
        add_k=convert(Uint64,LCG_ADD64)
        ran=convert(Uint64,1)
        un=convert(Uint64,n)
        while un>0
        	if un&1>0
                   ran=mul_k*ran+add_k
                end
        	add_k*=mul_k+1
               	mul_k*=mul_k
		un>>>=1
        end
        ran
end

@everywhere function updateArray(D,A::Array, count::Int64)
	global logTableSize;
	localTableSize=length(localpart(D))
	logLocalTableSize=convert(Uint64,log2(localTableSize))
	startGlobalAddress=convert(Uint64,(myid()-2))<<logLocalTableSize
	for i in 1:count
		localOffset=(A[i]>>>convert(Uint64,(64-logTableSize)))-startGlobalAddress+1
		(localpart(D))[localOffset]$=A[i]
	end
end


 @everywhere function HPCC_Power2NodesRandomAccess(D)
               sendCount=0
               global maxPendingUpdates+=0
               global logTableSize+=0
               global LCG_MUL64+=0
               global LCG_ADD64+=0
               global numProcs+=0
               totalPending=0
               pq=Collections.PriorityQueue()
               global Buckets+=0
               myID=convert(Uint64,myid()-2)   # 0-based
               localTableSize=length(localpart(D))
               logLocalTableSize=convert(Uint64,log2(localTableSize))
               startGlobalAddress=myID<<logLocalTableSize
               Ran=convert(Uint64,HPCC_starts(4*startGlobalAddress))
               while sendCount<(4*localTableSize)
                       if totalPending==maxPendingUpdates
                               X=Collections.peek(pq)
                               pq[X[1]]=0
                               ref=@spawnat X[1] updateArray(D,Buckets[X[1]-1,:],-1*X[2])
			       totalPending+=X[2]
			       #Buckets[X[1]-1,:]=zeros(maxPendingUpdates)-1
                       else
                               Ran=LCG_MUL64*Ran+LCG_ADD64
                               whichPE=(Ran>>>convert(Uint64,(64-logTableSize+logLocalTableSize)))&convert(Uint64,(numProcs-1))
                               if(whichPE==myID)
                                       localOffset=(Ran>>>convert(Uint64,(64-logTableSize)))-startGlobalAddress+1
                                       (localpart(D))[localOffset]$=Ran
                               else
                                       if Collections.haskey(pq,whichPE+2) 
						X=-1*pq[whichPE+2]+1 
				       else 
						pq[whichPE+2]=0
						X=1 
				       end
                                       pq[whichPE+2]-=1
                                       Buckets[whichPE+1,X]=Ran
                                       totalPending+=1
                               end
                               sendCount+=1
                       end
               end
	       while ~isempty(pq) && (Collections.peek(pq))[2]!=0
       			X=Collections.peek(pq)
       			ref=@spawnat X[1] updateArray(D,Buckets[X[1]-1,:],-1*X[2])
                        wait(ref)
			Collections.dequeue!(pq)
       	       end	
       end


   # Perform updates to main table.  The scalar equivalent is:
   #
   #     u64Int Ran;
   #     Ran = 1;
   #     for (i=0; i<NUPDATE; i++) {
   #       Ran = LCG_MUL64 * Ran + LCG_ADD64;
   #       Table[Ran >> (64 - LOG2_TABSIZE)] ^= Ran;
   #     }
   #


function verify(D)
	global TableSize;
	global LCG_MUL64;
	global LCG_ADD64;
	Table=Array(Uint64,TableSize)
	for i in 1:TableSize
		Table[i]=i-1
	end
	Ran=convert(Uint64,1)
	for i in 1:(TableSize*4)
		Ran=LCG_MUL64*Ran+LCG_ADD64
		Table[(Ran>>>convert(Uint64,64-log2(TableSize)))+1]$=Ran
	end
	
	count=0
	for i in 1:TableSize
		if D[i]!=Table[i]
			count+=1
		end
	end
	return count
end

A=@elapsed @sync {(@spawnat p HPCC_Power2NodesRandomAccess(HPCC_Table)) for p in procs(HPCC_Table)};
println("Elapsed time=",A)
Perf=(TableSize*4)/(1000000000*A)
println("Performance=",Perf," Gups")
HPCC_Table=convert(Array,HPCC_Table)
B=verify(HPCC_Table)
println("Failures=",B)

