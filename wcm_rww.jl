# Pkg.activate(".")

using OrdinaryDiffEq
using Plots

r_max = 500
aₑ = 310 #310e-9
bₑ = 125
dₑ = 0.16

aᵢ = 615 #615e-9
bᵢ = 177
dᵢ = 0.087

τₑ = 0.1
τᵢ = 0.01
γₑ = 0.641
γᵢ = 1
# const Iᵢ = 0.1e-9

function Hₑ(x)
    num_Hₑ = r_max + ((aₑ*x - bₑ - r_max)/(1 - exp(dₑ * (aₑ*x - bₑ - r_max))))
    den_Hₑ = 1 - exp(-dₑ * (aₑ*x - bₑ))
    num_Hₑ/den_Hₑ
end

function Hᵢ(x)
    num_Hᵢ = r_max + ((aᵢ*x - bᵢ - r_max)/(1 - exp(dᵢ * (aᵢ*x - bᵢ - r_max))))
    den_Hᵢ = 1 - exp(-dᵢ * (aᵢ*x - bᵢ))
    num_Hᵢ/den_Hᵢ
end

function wcm_rww_local_func!(du, u, p, t)
    Sₑ,Sᵢ= u
    ωₑₑ,ωₑᵢ,ωᵢₑ,ωᵢᵢ,Iₑ,Iᵢ = p
    @inbounds begin
        du[1] = dSₑ = (- Sₑ / τₑ) + (1 - Sₑ) * γₑ * Hₑ(ωₑₑ*Sₑ - ωᵢₑ*Sᵢ + Iₑ)
        du[2] = dSᵢ = (- Sᵢ / τᵢ) + (1 - Sᵢ) * γᵢ * Hᵢ(ωₑᵢ*Sₑ - ωᵢᵢ*Sᵢ + Iᵢ)
    end
    nothing
end

u0 = [0.7, 0.9] # Realistic range for Sₑ,Sᵢ - [0,1]
tspan = (0.0,100.0)
p = [0.21, 0.15, 1, 1, 0.382, 0.267]
prob = ODEProblem(wcm_rww_local_func!,u0,tspan,p)
sol = solve(prob, Tsit5())


pl1 = plot(sol[2000:3000], vars=(0,1), plotdensity=10000,lw=1.5)
pl2 = plot(sol[2000:3000], vars=(0,2), plotdensity=10000,lw=1.5)
pl3 = plot(sol[2000:3000], plotdensity=10000,lw=1.5)
