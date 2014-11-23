#!/nethome/kagarwal39/julia-0.3.1/julia/julia

function distribute(a::AbstractArray, pids, dist)
  owner = myid()
  rr = RemoteRef()
  put!(rr, a)
  DArray(size(a), pids, dist) do I
    remotecall_fetch(owner, ()->fetch(rr)[I...])
  end
end

ROW=1000;
COL=1000;
a=rand(ROW,COL);
b=Array(Float64,(ROW,COL));
t=Array(Float64,(ROW,COL));
ans=Array(Float64,(ROW,COL));
p=5;
q=5;
addprocs(p*q);
r=ifloor(ROW/p);
c=ifloor(COL/q);
for i=1:ROW
  for j=1:COL
    px=mod(i-1,p);
    py=mod(j-1,q);
    bx=ifloor((i-1)/p);
    by=ifloor((j-1)/q);
    t[px*r+bx+1,py*c+by+1]=a[i,j];
  end
end

#println(a);
a=t;
a=distribute(a,[2:nworkers()+1],[p,q]);
b=distribute(b,[2:nworkers()+1],[p,q]);

@everywhere function store(b::Array,a::Float64,x::Int64,y::Int64)
  b[x+1,y+1]=a; 
end

@everywhere function send(a::Array,b::DArray,p::Int64,q::Int64)
  id=myid();
  py=ifloor((id-2)/p); #minus 2 because worker ids begin from 2
  px=mod(id-2,p);
  for i=1:size(a)[1]
    for j=1:size(a)[2]

      #calculate new ids for transpose
      new_y=(i-1)*p+px+1;
      new_x=(j-1)*q+py+1;

      #new worker ids
      new_px=mod(new_x-1,p);
      new_py=mod(new_y-1,q);

      #new indexes
      new_bx=ifloor((new_x-1)/p);
      new_by=ifloor((new_y-1)/q);

      next_pid=new_px+new_py*p+2;
      #println("id=",myid()," next=",next_pid," new_bx=",new_bx," new_by=",new_by);
      @spawnat next_pid store(localpart(b),a[i,j],new_bx,new_by);
    end
  end
end

@time @sync { (@spawnat proc send(localpart(a),b,p,q)) for proc=procs(a) }

#println(b);

for i=1:ROW
  for j=1:COL
    px=mod(i-1,p);
    py=mod(j-1,q);
    bx=ifloor((i-1)/p);
    by=ifloor((j-1)/q);
    ans[i,j]=b[px*r+bx+1,py*c+by+1];
  end
end

#println(ans)
