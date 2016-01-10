#borrow from julia class

function mkmat_incidence_factor(b)
    factor=unique(b)
    coMat= spzeros(length(b),length(factor))

    dictFactor = Dict()
    index=1
    for i in factor
        dictFactor[i]=index
        index+=1
    end

    nrow=1
    for i in b
        myindex=dictFactor[i]
        coMat[nrow,myindex]=1
        nrow=nrow+1
    end
    return full(coMat),factor
end
