using ElectronGas
using GreenFunc
using Lehmann
using CompositeGrids
using Printf

rs = 5.0
beta = 4000.0
mass2 = 0.01

para = Parameter.rydbergUnit(1 / beta, rs, 3, Λs=mass2)

println(para)

function density(g)
    g = toTau(g)
    gk = g.dynamic[1, 1, :, end]
    for (ki, k) in enumerate(g.spaceGrid)
        gk[ki] = 4π * k^2 * dynamic(g, para.β - 1e-8, k, 1, 1)
    end
    n = Interp.integrate1D(gk, sigma.spaceGrid) / (2π)^3 * 2
    println("density = ", n)
    return n
end

sigma = SelfEnergy.G0W0(para, Nk=16, maxK=6 * para.kF)

println("mu: ", SelfEnergy.chemicalpotential(para, sigma))

sigma = toMatFreq(sigma)
mui = instant(sigma, para.kF, 1, 1)
println("static mu: ", mui)
mud = dynamic(sigma, 1, para.kF, 1, 1)
println("dynamic mu: ", mud)
ngrid = sigma.timeGrid
kgrid = sigma.spaceGrid

# ds_dw = dynamic(sigma, 1, para.kF, 1, 1)
# ds_dw = sigma.dynamic[1, 1, :,]
kF_label = searchsortedfirst(kgrid.grid, para.kF)
Σ_freq = GreenFunc.toMatFreq(sigma, [0, 1])

ΣI = imag(Σ_freq.dynamic[1, 1, kF_label, :])
Z0 = 1 / (1 - (ΣI[2] - ΣI[1]) / 2 / π * para.β)
println("z = ", Z0)

ds_dw = (Σ_freq.dynamic[1, 1, :, 2] - Σ_freq.dynamic[1, 1, :, 1]) / 2 / π * para.β
for (ki, k) in enumerate(kgrid)
    Ek = k^2 / (2 * para.me) - para.EF
    ds_dw[ki] *= 4π * k^2 * Spectral.fermiDirac(Ek, para.β)
end

# println()

println("∫dk dΣ/dω*f(ω) = ", Interp.integrate1D(ds_dw, kgrid) / (2π)^3 * 2)


g1 = deepcopy(sigma)
printstyled("g\n", color=:yellow)
for (ni, n) in enumerate(ngrid)
    for (ki, k) in enumerate(kgrid)
        Ek = k^2 / (2 * para.me) - para.EF
        g = Spectral.kernelFermiΩ(n, Ek, para.β)
        g1.dynamic[1, 1, ki, ni] = g
    end
end
n1 = density(g1)
println()

gg = deepcopy(sigma)
printstyled("g*g\n", color=:yellow)
for (ni, n) in enumerate(ngrid)
    for (ki, k) in enumerate(kgrid)
        Ek = k^2 / (2 * para.me) - para.EF
        g = Spectral.kernelFermiΩ(n, Ek, para.β)
        gg.dynamic[1, 1, ki, ni] = g * g
    end
end
n2 = density(gg)
println()

ggg = deepcopy(sigma)
printstyled("g*g*g\n", color=:yellow)
for (ni, n) in enumerate(ngrid)
    for (ki, k) in enumerate(kgrid)
        Ek = k^2 / (2 * para.me) - para.EF
        g = Spectral.kernelFermiΩ(n, Ek, para.β)
        ggg.dynamic[1, 1, ki, ni] = g * g * g
    end
end
n3 = density(ggg)
println()

gsg = deepcopy(sigma)
printstyled("g*sigma*g with static sigma\n", color=:yellow)
for (ni, n) in enumerate(ngrid)
    for (ki, k) in enumerate(kgrid)
        Ek = k^2 / (2 * para.me) - para.EF
        g = Spectral.kernelFermiΩ(n, Ek, para.β)
        # gsg.dynamic[1, 1, ki, ni] = g * g * sigma.dynamic[1, 1, ki, ni]
        gsg.dynamic[1, 1, ki, ni] = g * g * sigma.instant[1, 1, ki]
    end
end
println("$(density(gsg)), expect: $(n2*mui)\n")
println()

gsgsg = deepcopy(sigma)
printstyled("g*sigma*g sigma g with static sigma\n", color=:yellow)
for (ni, n) in enumerate(ngrid)
    for (ki, k) in enumerate(kgrid)
        Ek = k^2 / (2 * para.me) - para.EF
        g = Spectral.kernelFermiΩ(n, Ek, para.β)
        gsgsg.dynamic[1, 1, ki, ni] = g * sigma.instant[1, 1, ki] * g * sigma.instant[1, 1, ki] * g
    end
end
printstyled("$(density(gsgsg)), expect: $(n3*mui^2)\n", color=:red)
println()


gsg = deepcopy(sigma)
printstyled("g*sigma*g with dynamic sigma\n", color=:yellow)
for (ni, n) in enumerate(ngrid)
    for (ki, k) in enumerate(kgrid)
        Ek = k^2 / (2 * para.me) - para.EF
        g = Spectral.kernelFermiΩ(n, Ek, para.β)
        gsg.dynamic[1, 1, ki, ni] = g * g * sigma.dynamic[1, 1, ki, ni]
    end
end
println("$(density(gsg)), expect: $(n2*mud)\n")
println()

gsgsg = deepcopy(sigma)
printstyled("g*sigma*g*sigma*g with dynamic sigma\n", color=:yellow)
for (ni, n) in enumerate(ngrid)
    for (ki, k) in enumerate(kgrid)
        Ek = k^2 / (2 * para.me) - para.EF
        g = Spectral.kernelFermiΩ(n, Ek, para.β)
        gsgsg.dynamic[1, 1, ki, ni] = g * sigma.dynamic[1, 1, ki, ni] * g * sigma.dynamic[1, 1, ki, ni] * g
    end
end
printstyled("$(density(gsgsg)), expect: $(n3*mud^2)\n", color=:red)
println()

# sigma = toTau(sigma)
# println(dynamic(sigma, para.β - 1e-6, para.kF, 1, 1))


