module QTL

using DataFrames

include("genotypes.jl")
include("fixed_effects.jl")
include("QTL_types.jl")
include("tools.jl")
include("files.jl")
include("qc.jl")
include("Gibbs.jl")

end
