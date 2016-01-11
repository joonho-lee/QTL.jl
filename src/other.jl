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
                                 marker_matrix::MarkerMatrix,
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
