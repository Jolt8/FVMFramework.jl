include(joinpath(@__DIR__, "FVMVerification.jl"))
using .FVMVerification

if isempty(ARGS)
    level = :unit
else
    level = Symbol(ARGS[1])
end

if length(ARGS) >= 2
    output_directory = ARGS[2]
else
    output_directory = joinpath(@__DIR__, "reports", string(level))
end

report = validate(; level = level, output_directory = output_directory)
if has_failures(report)
    exit(1)
end
