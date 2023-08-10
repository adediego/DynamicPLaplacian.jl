using DifferentialEquations 
using Memoization
using LinearAlgebra
using ForwardDiff

function rot_double_gyre!(du, u, p, t)
    du1_P =   π * sin(2π*u[1]) * cos(π*u[2])
    du2_P = -2π * cos(2π*u[1]) * sin(π*u[2])

    du1_F =  2π * sin(π*u[1]) * cos(2π*u[2])
    du2_F =  -π * cos(π*u[1]) * sin(2π*u[2])

    s = t^2 * (3-2t)
    if t<0
        s = zero(s)
    elseif t>1
        s = one(s)
    end

    du[1] = (1-s) * du1_P + s * du1_F
    du[2] = (1-s) * du2_P + s * du2_F
end

@memoize Dict function rot_gyre_flow(u0, T)
    tspan = (0.0, T)
    u0_array = [u0[1], u0[2]]
    prob = ODEProblem(rot_double_gyre!, u0_array, tspan)
    sol_array = DifferentialEquations.solve(prob, Tsit5(), abstol=1e-7, reltol=1e-7)(T)
    return VectorValue(sol_array[1], sol_array[2])
end

@memoize Dict function dinv_flow(u0, flow; diffmethod = :finite_diff)

    T = 1.0
    if diffmethod == :finite_diff
        δ = 1e-4
        du1 = VectorValue(δ, 0.0)
        du2 = VectorValue(0.0, δ)

        dflow1 = (flow(u0 .+ du1/2, T) .- flow(u0 .- du1/2, T)) / δ
        dflow2 = (flow(u0 .+ du2/2, T) .- flow(u0 .- du2/2, T)) / δ

        D = TensorValue{2,2,eltype(u0)}(dflow1[1], dflow1[2], dflow2[1], dflow2[2])
        Dinv = inv(D)
        return Dinv
    elseif diffmethod == :forwarddiff
        flow_func = x -> [flow(x, T)...]
        jac = ForwardDiff.jacobian(flow_func, [u0...])
        D = TensorValue{2,2,eltype(u0)}(jac[1,1], jac[2,1], jac[1,2], jac[2,2])
        return inv(D)
    else
        throw("Unknown diffmethod")
    end
 end

function CG_tensor(u0, flow; diffmethod=:finite_diff)
    Dinv = dinv_flow(u0, flow; diffmethod)
    return Dinv ⋅ Dinv'
end

@memoize Dict function cg_eigmax(u0, flow; diffmethod=:finite_diff)
    dd = CG_tensor(u0, flow; diffmethod)
    A = [dd[1,1] dd[1,2]; dd[2,1] dd[2,2]]
    return eigmax(A)
end

@memoize Dict function mean_CG_tensor(u0, flow; diffmethod=:finite_diff)
    Cinv = CG_tensor(u0, flow; diffmethod)
    id = TensorValue(1.0, 0.0,0.0, 1.0)
    return 0.5(id + Cinv)
end



