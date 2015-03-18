function read_genotypes(file,nrow,ncol,header=true)
    f=open(file)
    mat = []

    if header==true
        readline(f)
        mat = zeros(Int64,nrow-1,ncol)
    else
        mat = zeros(Int64,nrow,ncol)
    end

    for i=1:(nrow-1)
        mat[i,:]=int(split(readline(f)))

        if(i%1000==0)
            println("This is line ",i)
        end
    end

    close(f)

    return mat
end


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

function findfixedloci(X)
    nrow,ncol=size(X)
    fixedloci=Int64[]

    for i = 1:ncol
        if var(X[:,i])==0.0
            push!(fixedloci,i)
        end

        if(i%10000==0)
            println("This is column ",i)
        end
    end

    return fixedloci
end



#function rw2binary()
