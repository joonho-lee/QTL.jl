type InputParameters2
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


type Output
    meanFxdEff::Array{Float64,1}
    meanMrkEff::Array{Float64,2}
    mdlFrq::Array{Float64,1}
    resVar::Array{Float64,1}
    genVar::Array{Float64,1}
    pi::Array{Float64,1}
    scale::Array{Float64,1}

    function Output(input::InputParameters,random::GibbsMats,fixed::GibbsMats)
      nFixedEffects = fixed.ncols
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


