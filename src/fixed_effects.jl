type FixedMatrix
    C::Array{Float64,2}
    variables::Array{Any,1}
end

function make_fixed(file;id4col=true,typeofvariable=true) #maybe better to use number, factor
    myfile = open(file)
    #get number of columns
    row1   = split(readline(myfile))

    if id4col==true
      fixedNames=row1 #skip header
    else
      fixedNames= NaN
    end
    #set types for each column and get number of markers
    ncol= length(row1)
    if typeofvariable==true
        etv = Array(DataType,ncol)
        fill!(etv,UTF8String)
    else
        etv = typeofvariable
    end
    close(myfile)

    #read genotypes
    df = readtable(file, eltypes=etv, separator = ' ',header=id4col)

    variables=Array{Any}(ncol)
    if typeofvariable[1]==Float64
        C = convert(Array,df[:,1])
        variables[1] = "covariate"
    else
        C,variables[1]=mkmat_incidence_factor(df[:,1])
    end
    for i=2:size(df,2)
      if typeofvariable[i]==Float64
          C = [C convert(Array,df[:,i])]
          variables[i] = "covariate"
      else
          C,variables[i]=[C mkmat_incidence_factor(df[:,i])]
      end
    end

    return FixedMatrix(C,variables)
end

function mkmat_incidence_factor(b)
    factor=unique(b)
    coMat= spzeros(length(b),length(factor)) #maybe better another way to construct sparse matrix here

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
    return coMat,factor
end

export FixedMatrix
#=
    ids = Array(UTF8String,size(fixed.C,2))
    k = 1
    for i in fixed.variable
      if variable[i]=="covariate"
        ids[k] = "covariate"
        k+=1
      else
        for level in variabel[i]
          ids[k]=level
=#

