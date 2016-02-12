type InputParameters
  method::AbstractString #BaysABC
  chainLength::Int64     # number of iterations
  probFixed::Float64     # parameter "pi" the probability SNP effect is zero
  varGenotypic::Float64  # used to derive hyper parameter (scale) for locus effect variance
  varResidual::Float64   # used to derive hyper parameter (scale) for locus effect variance
  estimateVariance::Bool # true or false
  estimatePi::Bool       # true or false
  estimateScale::Bool    # true or false
  dfEffectVar::Float64   # hyper parameter (degrees of freedom) for locus effect variance
  nuRes::Float64         # hyper parameter (degrees of freedom) for residual variance
  nuGen::Float64         # hyper parameter (degree of freedom) for genetic variance (for ϵ)
end

InputParameters()=InputParameters("BayesC",50000,0.95,1.0,1.0,true,false,false,4,4,4)

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
    varGenotypic::Float64      #genotypic variance
    varResidual::Float64       #residual variance
    varEffect::Float64  # common marker variance
    scaleVar::Float64   # scale factor for locus effects
    scaleRes::Float64   # scale factor for residual varianc
    scaleGen::Float64
    fixed_effects       # sample of fixed effects β
    α                   # sample of partial marker effects unconditional on δ
    δ                   # inclusion indicator for marker effects
    u                   # sample of marker effects,u=α·δ
    π                   # probFixed
    locusEffectVar      #locus-specific variance

    iter
    yCorr

    imputation_residual #residual in SSBR

    function Current(input::InputParameters,geno::Genotypes,
                     fixed::FixedMatrix,y::Array{Float64,1})
        varGenotypic= input.varGenotypic
        varResidual = input.varResidual
        π           = input.probFixed
        dfEffectVar = input.dfEffectVar
        nuRes       = input.nuRes
        nuGen       = input.nuGen
        varGenotypic= input.varGenotypic
        sum2pq      = geno.sum2pq
        nMarkers    = geno.nMarkers
        nFixedEffects  = size(fixed.C,2)

        varEffect      = varGenotypic/((1-π)*sum2pq)
        locusEffectVar = fill(varEffect,nMarkers) #add if statement later
        scaleVar       = varEffect*(dfEffectVar-2)/dfEffectVar # scale factor for locus effects
        scaleRes       = varResidual*(nuRes-2)/nuRes                   # scale factor for residual varianc
        scaleGen       = varGenotypic*(nuGen-2)/nuGen                   # scale factor for residual varianc
        β          = zeros(nFixedEffects)  # sample of fixed effects
        α          = zeros(nMarkers)       # sample of partial marker effects unconditional on δ
        δ          = zeros(nMarkers)       # inclusion indicator for marker effects
        u          = zeros(nMarkers)       # sample of marker effects

        iter       = 0
        yCorr      = copy(y)
        ϵ          = zeros(1)

        new(varGenotypic,varResidual,varEffect,scaleVar,scaleRes,scaleGen,
            β,α,δ,u,π,locusEffectVar,iter,yCorr,ϵ)
    end
end

type Output
    mean_fixed_effects::Array{Float64,1}
    meanMarkerEffects::Array{Float64,1}
    modelFreq::Array{Float64,1}
    resVar::Array{Float64,1}
    genVar::Array{Float64,1}
    pi::Array{Float64,1}
    scale::Array{Float64,1}

    mean_imputation_residual::Array{Float64,1}

    function Output(input::InputParameters,geno::Genotypes,fixed::FixedMatrix)
      chainLength   = input.chainLength
      nFixedEffects = length(fixed.variables)
      nMarkers      = geno.nMarkers
      new(zeros(nFixedEffects),
          zeros(nMarkers),
          zeros(nMarkers),
          zeros(chainLength),
          zeros(chainLength),
          zeros(nMarkers),
          zeros(chainLength),
          zeros(1))
    end
end

export InputParameters
export Output
export GibbsMats
export Current
