include(joinpath(@__DIR__, "slau2", "SLAU2Verification.jl"))
using .SLAU2Verification

if isempty(ARGS)
    level = :unit
else
    level = Symbol(ARGS[1])
end

if length(ARGS) >= 2
    output_directory = ARGS[2]
else
    output_directory = joinpath(@__DIR__, "reports", "slau2", string(level))
end

report = validate_slau2(; level = level, output_directory = output_directory)
if has_failures(report)
    exit(1)
end
