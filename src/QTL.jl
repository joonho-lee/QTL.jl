module QTL

using DataFrames
using Distributions

include("genotypes.jl")
include("fixed_effects.jl")
include("QTL_types.jl")
include("tools.jl")
include("files.jl")
include("QualityControl.jl")
include("Gibbs.jl")

end
