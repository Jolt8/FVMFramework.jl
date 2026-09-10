using Documenter
using FVMFramework # Replace with your actual package name

makedocs(
    sitename = "FVMFramework",
    modules = [FVMFramework],
    authors = "Will Martin <willemartin@icloud.com>",
    pages = [
        "FVMFramework.jl: A Framework for solving Finite Volume Method Problems in Julia" => "index.md",
        "Setting up a Simulation" => "quickstart/heat_transfer_example.md",
        "Useful Features" => "useful_features/for_fields_and_foreach_field_at.md",
        "Included Flux Physics" => "included_physics/flux_physics/advection.md",
        "Included Flux Physics" => "included_physics/flux_physics/darcy_flow.md",
        "Included Flux Physics" => "included_physics/flux_physics/diffusion.md",
        "Included Physics" => "included_physics/flux_physics/heat_transfer.md"
    ]
)

deploydocs(
    repo = "https://github.com/Jolt8/FVMFramework.jl",
)