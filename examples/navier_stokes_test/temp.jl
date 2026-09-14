
test = ComponentVector(
    zoop = 1,
)

test = merge_properties(test, ComponentVector(zap = [1.0, 1.0, 1.0]))

test.zap

test_2 = ComponentVector(
    velocity = [
        u_named[1].u_vel, 
        u_named[1].v_vel, 
        u_named[1].w_vel
    ]
)

velocity = []

for i in eachindex(u_named[1].u_vel)
    push!(velocity, collect([u_named[1].u_vel[i], u_named[1].v_vel[i], u_named[1].w_vel[i]]))
end

velocity

test_2 = ComponentVector(velocity = velocity,)

mat = reduce(hcat, velocity)

mat = transpose(mat)

test = merge_properties(test, ComponentVector(velocity = mat))
