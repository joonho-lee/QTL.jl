#!/usr/bin/env julia

function missing2mean(X,missing=9)
    nrow,ncol = size(X)
    for i=1:ncol
       index=find(x->x==missing,X[:,i])
       cols = [1:nrow]
       deleteat!(cols,index)
       X[index,i]=int(mean(X[cols,i]))

       if(i%3000==0)
          println("This is line ",i)
       end
     end
end


file=readdlm(ARGS[1])
missing2mean(file)
writedlm(ARGS[2],file)

