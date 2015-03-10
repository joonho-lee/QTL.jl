module QTL


function get_column(X,j)
  nrow,ncol = size(X)
  if j>ncol||j<0
      error("column number is wrong!")
  end
  indx = 1 + (j-1)*nrow
  ptr = pointer(X,indx)
  pointer_to_array(ptr,nrow)
end

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
        myindex=DictFactor[i]
        coMat[nrow,myindex]=1
        nrow=nrow+1
    end
    return full(coMat),dictFactor
end

function delete_rc(x,rowi,colj)
    numRow,numCol = size(x)
    nrow = [1:numRow]
    ncol = [1:numCol]

    if j != 0
        deleteat!(nrow,rowi)
    end

    if k != 0
        deleteat!(ncol,colj)
    end

    return x[nrow,ncol]
end

function find_column_byname(file,names)
  colname = split(readline(file))
  nameIndex=Int64[]

  for i in names
    index=find(colname,names)
    push!(nameIndex,index)
  end

  return nameIndex
end

end