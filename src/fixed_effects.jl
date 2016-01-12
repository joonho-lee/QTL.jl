type FixedMatrix
    C::Array{Float64,2}
    nFixedEffects = size(C,2)
    fixedNames::Array{AbstractString,1}

    function FixedMatrix(file::DataFrames) ###More
        mean2pq     = (2*p*(1-p)')[1,1]
        xArray = get_column_ref(X)
        XpRinvX = getXpRinvX(X)
        new(X,xArray,XpRinvX,markerMeans,mean2pq,true)
    end
end

function make_fixed(file;id4col=true,variable=true)
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
    if variable==true
        etv = Array(DataType,ncol)
        fill!(etv,UTF8String)
    else
        etv = variable
    end
    close(myfile)

    #read genotypes
    df = readtable(file, eltypes=etv, separator = ' ',header=id4col)

    if variable[1]==Float64
        C = convert(Array,df[:,1])
    else
        C,levels=mkmat_incidence_factor(df[:,1])
    end
    for i=2:size(df,2)
      if variable[i]==Float64
          C = [C convert(Array,df[:,i])]
      else
          C,levels=[C mkmat_incidence_factor(df[:,i])]
      end
    end

    return
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
        myindex=dictFactor[i]
        coMat[nrow,myindex]=1
        nrow=nrow+1
    end
    return full(coMat),factor
end
