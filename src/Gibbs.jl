function sample_fixed!(mats::GibbsMats,current::Current,out::Out)
    nObs,nEffects = mats.nrows,mats.ncols
    xArray        = mats.XArray
    xpx           = mats.xpx
    yCorr         = current.yCorr
    α             = current.α
    varRes        = current.varResidual
    iIter         = 1/current.iter
    meanAlpha     = out.meanAlpha

    for j=1:nEffects
        x = xArray[j]
        rhs = dot(x,yCorr) + xpx[j]*α[j,1]
        lhs      = xpx[j]
        invLhs   = 1.0/lhs
        mean     = invLhs*rhs
        oldAlpha = α[j,1]
        α[j]     = mean + randn()*sqrt(invLhs*varRes)
        BLAS.axpy!(oldAlpha-α[j,1],x,yCorr)
        meanAlpha[j] += (α[j] - meanAlpha[j])*iIter
    end
end

function sample_random_ycorr!(mats::GibbsMats,current::Current,out::Out)#sample vare and vara
    nObs,nEffects = mats.nrows,mats.ncols
    xArray        = mats.XArray
    xpx           = mats.xpx
    yCorr         = current.yCorr
    α             = current.α
    varRes        = current.varResidual
    λ             = current.varResidual/current.varEffects
    iIter         = 1/current.iter
    meanAlpha     = out.meanAlpha

    for j=1:nEffects
        x = xArray[j]
        rhs = dot(x,yCorr) + xpx[j]*α[j,1]
        lhs      = xpx[j] + λ
        invLhs   = 1.0/lhs
        mean     = invLhs*rhs
        oldAlpha = α[j,1]
        α[j]     = mean + randn()*sqrt(invLhs*varRes)
        BLAS.axpy!(oldAlpha-α[j,1],x,yCorr)
        meanAlpha[j] += (α[j] - meanAlpha[j])*iIter
    end
end

function sample_random_rhs!(lhs,rhs,current::Current,out::Out) #(Gianola Book)
    nEffects = size(lhs,1)
    α             = current.α
    varRes        = current.varResidual
    iIter         = 1/current.iter
    meanAlpha     = out.meanAlpha

    for (i in 1:nEffects) #argument lhs here is a sparse matrix for whole lhs
        α[i] = 0.0
        rhsi = rhs[i] - lhs[i,:]*α
        lhsi = lhs[i,i]
        invLhs = 1.0/lhsi
        meani  = invLhs*rhsi[1]
        α[i] = meani + randn()*sqrt(invLhs*varRes)
        meanAlpha[i] += (α[i] - meanAlpha[i])*iIter
    end
end


## NOT sample vare and vara
## construct lhsDi,sd outside

function sample_random_ycorr!(mats::GibbsMats,current::Current,out::Out,lhsDi,sd)#not sample vare and vara
    nObs,nEffects = mats.nrows,mats.ncols
    xArray        = mats.XArray
    xpx           = mats.xpx
    yCorr         = current.yCorr
    α             = current.α
    iIter         = 1/current.iter
    meanAlpha     = out.meanAlpha

    for j=1:nEffects
        @inbounds x = xArray[j]
        @inbounds rhs = dot(x,yCorr) + xpx[j]*α[j,1]
        @inbounds mean     = lhsDi[j]*rhs
        @inbounds oldAlpha = α[j,1]
        @inbounds α[j]     = mean + randn()*sd[j]
        @inbounds BLAS.axpy!(oldAlpha-α[j,1],x,yCorr)
        @inbounds meanAlpha[j] += (α[j] - meanAlpha[j])*iIter
    end
end

function sample_random_rhs!(lhsCol,rhs,current::Current,out::Out,lhsDi,sd)
    nEffects = size(lhs,1)      #arguments lhs here is a array of cols of lhs
    α             = current.α
    iIter         = 1/current.iter
    meanAlpha     = out.meanAlpha

    for (i in 1:n)
        @inbounds b[i] = 0.0
        @inbounds rhsi = rhs[i] - lhsCol[i,:]'b
        @inbounds α[i] = lhsDi[i]*rhsi[1] + randn()*sd[i]
        @inbounds meanAlpha[i] += (α[i] - meanAlpha[i])*iIter
    end
end


type OutputValues
    meanFxdEff::Array{Float64,1}
    meanMrkEff::Array{Float64,2}
    mdlFrq::Array{Float64,1}
    resVar::Array{Float64,1}
    genVar::Array{Float64,1}
    pi::Array{Float64,1}
    scale::Array{Float64,1}

    function OutputValues(input_parameters,marker_matrix,fixed_matrix)
      nFixedEffects = fixed_matrix.nFixedEffects
      nMarkers      = marker_matrix.nMarkers
      chainLength   = input_parameters.chainLength
      new(zeros(nFixedEffects),
          zeros(nMarkers,1),
          zeros(nMarkers),
          zeros(chainLength),
          zeros(chainLength),
          zeros(chainLength),
          zeros(chainLength))
    end
end

type EstimatedParameters
    vare::Float64
    varEffects::Float64
    scaleVar::Float64   # scale factor for locus effects
    scaleRes::Float64   # scale factor for residual varianc
    β                   # sample of fixed effects
    α                   # sample of partial marker effects unconditional on δ
    u                   # sample of marker effects
    δ                   # inclusion indicator for marker effects
    π                   # probFixed
    locusEffectVar      #locuda-specific variance

    function EstimatedParameters(input_parameters::InputParameters,
                                 marker_matrix::MarkerMatrix
                                 fixed_matrix::FixedMatrix)
        vare        = input_parameters.varResidual
        π           = input_parameters.probFixed
        dfEffectVar = input_parameters.dfEffectVar
        nuRes       = input_parameters.nuRes
        varGenotypic= input_parameters.varGenotypic
        mean2pq     = marker_matrix.mean2pq
        nMarkers    = marker_matrix.nMarkers
        nFixedEffects  = fixed_matrix.nFixedEffects

        varEffects     = varGenotypic/((1-π)*mean2pq)
        locusEffectVar = fill(varEffects,nMarkers)
        scaleVar       = varEffects*(dfEffectVar-2)/dfEffectVar # scale factor for locus effects
        scaleRes       = vare*(nuRes-2)/nuRes                   # scale factor for residual varianc
        β          = zeros(nFixedEffects)  # sample of fixed effects
        α          = zeros(nMarkers)       # sample of partial marker effects unconditional on δ
        δ          = zeros(nMarkers)       # inclusion indicator for marker effects
        u          = zeros(nMarkers)       # sample of marker effects
        new(vare,varEffects,scaleVar,scaleRes,β,α,u,δ,π,locusEffectVar)
    end
end
