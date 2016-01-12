function sample_fixed!(mats::GibbsMats,current::Current,out::Output)
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

function sample_random_ycorr!(mats::GibbsMats,current::Current,out::Output)#sample vare and vara
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

function sample_random_rhs!(lhs,rhs,current::Current,out::Output) #(Gianola Book)
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

function sample_random_ycorr!(mats::GibbsMats,current::Current,out::Output,lhsDi,sd)#not sample vare and vara
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

function sample_random_rhs!(lhsCol,rhs,current::Current,out::Output,lhsDi,sd)
    nEffects = size(lhs,1)      #arguments lhs here is a array of cols of lhs
    α             = current.α
    iIter         = 1/current.iter
    meanAlpha     = out.meanAlpha

    for (i in 1:nEffects)
        @inbounds α[i] = 0.0
        @inbounds rhsi = rhs[i] - lhsCol[i]'α
        @inbounds α[i] = lhsDi[i]*rhsi[1] + randn()*sd[i]
        @inbounds meanAlpha[i] += (α[i] - meanAlpha[i])*iIter
    end
end
