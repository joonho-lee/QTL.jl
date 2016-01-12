type Genotypes
  obsID::Array{UTF8String,1}  #row ID of genotypes
  markerID
  nObs::Int64
  nMarkers::Int64
  alleleFreq::Array{Float64,2}
  sum2pq::Float64
  centered::Bool
  genotypes::Array{Float64,2}
end

function make_genotypes(file;id4col=false,center=true)
    myfile = open(file)
    #get number of columns
    row1   = split(readline(myfile))

    if id4col==true
      markerID=row1[2:end]  #skip header
    else
      markerID= NaN
    end
    #set types for each column and get number of markers
    ncol= length(row1)
    etv = Array(DataType,ncol)
    fill!(etv,Float64)
    etv[1]=UTF8String
    close(myfile)

    #read genotypes
    df = readtable(file, eltypes=etv, separator = ' ',header=id4col)

    obsID=convert(Array,df[:,1])
    genotypes =convert(Array,df[:,2:end])
    nObs,nMarkers = size(genotypes)

    if center==true
        markerMeans = center!(genotypes) #centering
    else
        markerMeans = center(genotypes)  #get marker means
    end
    p             = markerMeans/2.0
    mean2pq       = (2*p*(1-p)')[1,1]

    return Genotypes(obsID,markerID,nObs,nMarkers,p,mean2pq,center,genotypes)
end

export make_genotypes
export Genotypes
