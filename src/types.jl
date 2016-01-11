type InputParameters
  method::AbstractString #BaysABC
  chainLength::Int64     # number of iterations
  probFixed::Float64     # parameter "pi" the probability SNP effect is zero
  varGenotypic::Float64  # used to derive hyper parameter (scale) for locus effect variance
  varResidual::Float64   # used to derive hyper parameter (scale) for locus effect variance
  estimateVariance::bool # true or false
  estimatePi::bool       # true or false
  estimateScale::bool    # true or false
  dfEffectVar::Float64   # hyper parameter (degrees of freedom) for locus effect variance
  nuRes::Float64         # hyper parameter (degrees of freedom) for residual variance
end

InputParameters()=InputParameters("BayesC",50000,0.95,1.0,1.0,true,false,false,4,4)


type GibbsMats
    X::Array{Float64,2}
    nrows::Int64
    ncols::Int64
    xArray::Array{Array{Float64,1},1}
    xpx::Array{Float64,1}
    function GibbsMats(X::Array{Float64,2}) ###More
        nrows,ncols = size(X)
        xArray = get_column_ref(X)
        XpX = getXpRinvX(X)
        new(X,nrows,ncols,xArray,XpX)
    end
end

type Current
  yCorr
  β
  α
  varEffects
  varResidual
  iter
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
    locusEffectVar      #locus-specific variance
    iter
    yCorr

    function EstimatedParameters(input::InputParameters,
                                 marker_matrix::MarkerMatrix
                                 fixed_matrix::FixedMatrix)
        vare        = input.varResidual
        π           = input.probFixed
        dfEffectVar = input.dfEffectVar
        nuRes       = input.nuRes
        varGenotypic= input.varGenotypic
        mean2pq     = marker_matrix.mean2pq
        nMarkers    = marker.ncols
        nFixedEffects  = fix.ncols

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


type Output
    meanFxdEff::Array{Float64,1}
    meanMrkEff::Array{Float64,2}
    mdlFrq::Array{Float64,1}
    resVar::Array{Float64,1}
    genVar::Array{Float64,1}
    pi::Array{Float64,1}
    scale::Array{Float64,1}

    function Output(input::InputParameters,random::GibbsMats,fix::GibbsMats)
      nFixedEffects = fix.ncols
      nMarkers      = random.ncols
      chainLength   = input.chainLength
      new(zeros(nFixedEffects),
          zeros(nMarkers,1),
          zeros(nMarkers),
          zeros(chainLength),
          zeros(chainLength),
          zeros(chainLength),
          zeros(chainLength))
    end
end


