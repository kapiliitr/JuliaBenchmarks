#!/nethome/kagarwal39/julia-0.3.1/julia/julia

const P=6;
const Q=16;
const ROW=2400;
const COL=2400;

isdefined(:P) || (const P = 2); #Number of process in x dimension
isdefined(:Q) || (const Q = 3); #Number of processes in y dimension
isdefined(:ROW) || (const ROW = 30); #Number of rows in matrix
isdefined(:COL) || (const COL = 30); #Number of columns in matrix
isdefined(:r) || (const r = 5); #Number of rows in a block
isdefined(:c) || (const c = 5); #Number of rows in a block

if r!=c
  println("Only square blocks supported yet.");
  exit(0);
end

if(mod(ROW,P*r)!=0 || mod(COL,Q*c)!=0)
  println("ERROR : Rows and columns are not exactly divisible into the blocks and workers.");
  exit(0);
end 

addprocs(P*Q);

#Customizing the distribution of matrix among the workers
function distribute(a::AbstractArray, pids, dist)
  owner = myid()
  rr = RemoteRef()
  put!(rr, a)
  DArray(size(a), pids, dist) do I
    remotecall_fetch(owner, ()->fetch(rr)[I...])
  end
end

function main()

  matA = rand(ROW,COL); #Initial matrix to be transposed
  matB = Array(Float64,(ROW,COL)); #Matrix containing the transpose
  matTemp = Array(Float64,(ROW,COL)); #Temproray matrix
  matT = Array(Float64,(ROW,COL));

  row_b = ifloor(ROW/r); #Number of blocks in x direction
  col_b = ifloor(COL/c); #Number of blocks in y direction
  b_x_t = ifloor(row_b/P); #Number of blocks per worker in x direction
  b_y_t = ifloor(col_b/Q); #Number of blocks per worker in y direction
  for i = 1:row_b
    for j = 1:col_b

      x_t = mod(i-1,P); #threadId.x
      y_t = mod(j-1,Q); #threadId.y
      x_b = ifloor((i-1)/P); #blockId.x
      y_b = ifloor((j-1)/Q); #blockId.y
      
      for m = 1:r
        for n = 1:c
          new_b_x = x_t * b_x_t + x_b;
          new_b_y =y_t * b_y_t + y_b;
          matTemp[new_b_x * r + m, new_b_y * c + n] = matA[(i-1) * r + m,(j-1) * c + n];
        end
      end

    end
  end

  matOrig = matA; #Saving original matrix
  matA = matTemp;
  matA = distribute(matA,[2:nworkers()+1],[P,Q]);
  matB = distribute(matB,[2:nworkers()+1],[P,Q]);

  @time @sync { (@spawnat proc dotrans(localpart(matA),matB,P,Q,r,c)) for proc = procs(matA) }

  for i = 1:row_b
    for j = 1:col_b

      x_t = mod(i-1,P);
      y_t = mod(j-1,Q);
      x_b = ifloor((i-1)/P);
      y_b = ifloor((j-1)/Q);
    
      for m = 1:c
        for n = 1:r
          new_b_x = x_t * b_x_t + x_b;
          new_b_y = y_t * b_y_t + y_b;
          matT[(i-1) * c + m,(j-1) * r + n] = matB[new_b_x * c + m, new_b_y * r + n];
        end
      end

    end
  end

  verifyTrans(matOrig,matT,ROW,COL);
 
end

function verifyTrans(matA::Array,matB::Array,R::Int64,C::Int64)
  correct = 1;  
  
  for i = 1:R
    for j = 1:C
      if matB[j,i] != matA[i,j]
        correct = 0;
        break;
      end
    end
    if correct == 0
      break;
    end
  end

  if correct == 1
    println("Result passed.");
  else
    println("Transpose failed.");
  end
end

@everywhere function store(matB::Array,trans::Array,b_x::Int64,b_y::Int64,r::Int64,c::Int64)
  for i = 1:r
    for j = 1:c
      matB[b_x * r + i,b_y * c + j] = trans[i,j]; 
    end
  end
end

@everywhere function dotrans(matA::Array,matB::DArray,P::Int64,Q::Int64,r::Int64,c::Int64)
  id = myid();
  y_t = ifloor((id-2)/P); #minus 2 because worker ids begin from 2
  x_t = mod(id-2,P);
  
  row_b = ifloor(size(matA)[1]/r); 
  col_b = ifloor(size(matA)[2]/c);

  for i = 1:row_b
    for j = 1:col_b

      #Locally transpose a block of matrix
      transTemp = Array(Float64,(c,r));
      for m = 1:r
        for n = 1:c
          transTemp[n,m] = matA[(i-1) * r + m,(j-1) * c + n];
        end
      end

      #calculate new ids for transpose
      new_y = (i-1) * P + x_t + 1;
      new_x = (j-1) * Q + y_t + 1;

      #new worker ids
      new_x_t = mod(new_x - 1,P);
      new_y_t = mod(new_y - 1,Q);

      #new indexes
      new_x_b = ifloor((new_x - 1)/P);
      new_y_b = ifloor((new_y - 1)/Q);
      
      next_pid = new_x_t + new_y_t * P + 2;
      @spawnat next_pid store(localpart(matB),transTemp,new_x_b,new_y_b,c,r);
    end
  end
end

# Calling the main function
main()
