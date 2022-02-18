# Calculate the susceptibility vs. matsubara frequency 
# using QuantumStatistics, LinearAlgebra, Printf, Parameters, DelimitedFiles
using Lehmann, LinearAlgebra, Printf, Parameters, DelimitedFiles
using Roots, Polylogarithms, FastGaussQuadrature
using Gaston, LaTeXStrings, ElectronGas
#using QuantumStatistics

include("parameter.jl")

function chemical_potential(beta)
    f(β, μ) = real(polylog(3 / 2, -exp(β * μ))) + 4 / 3 / π^0.5 * (β)^(3 / 2)
    g(μ) = f(beta, μ)
    return find_zero(g, (-10000, 1))
end

@with_kw struct Para
    Qsize::Int = 5
    nsize::Int = 11
    qGrid::Vector{Float64} = [q for q in LinRange(0.0, 0.1 * kF, Qsize)]
    nGrid::Vector{Int64} = [n for n in 0:nsize-1]
    μ = chemical_potential(beta) * EF
    χdlr = DLRGrid(EF * 100, β, 1e-8, false, :ph)
end

function susceptibility(para, QBin)
    #χ_tau = zeros(Float64, para.χdlr.size)
    χ0 = zeros(Float64, para.nsize)
    # χ = similar(χ0)
    q = para.qGrid[QBin]
    @printf("ExtMom: %f, DLR nsize: %d\n", q, para.χdlr.size)

    filename = "Data_Polar/polar_beta$(beta)_rs$(rs)_lam$(λ)_o$(order)_q$(QBin-1).dat1"
    print("Loading ", filename, "\n")
    polar_tau = readdlm(filename)
    χ_tau::Vector{Float64} = polar_tau[6:end, 2]
    #print(χ_tau)    
    #print(χ)
    χ = tau2matfreq(para.χdlr, χ_tau, para.nGrid)
    #χ = DLR.tau2matfreq(:corr, χ_tau, para.χdlr, para.nGrid, axis=1)
    #print(tau2matfreq(para.χdlr,χ_tau))

    param = Parameter.rydbergUnit(1.0 / (β * EF), rs, dim; μ = para.μ)

    # print(" n χ_0 χ for q/kF=$(q/param.kF)\n")
    # println(Polarization.Polarization0_FiniteTemp(q, 0, param) * 2)
    # println(Polarization.Polarization0_ZeroTemp(q, 0, param) * 2)
    for (idn, n) in enumerate(para.nGrid)
        # χ0[idn] = TwoPoint.LindhardΩnFiniteTemperature(param.dim, q, n, para.μ, param.kF, param.β, param.me, param.spin)[1]
        χ0[idn] = Polarization.Polarization0_FiniteTemp(q, n, param) * param.spin
        # @printf("%d %10.8f %10.8f \n", n, -χ0[idn], χ[idn])
    end
    return -χ0, real(χ), param
end

function run()
    para = Para()
    @unpack Qsize, qGrid, nGrid = para
    # p = Gaston.Figure[]
    wnGrid = nGrid * 2π / β
    for qi in 1:Qsize
        χ0, χ, param = susceptibility(para, qi)
        # y0 = χ0.*(wnGrid.*wnGrid)/qGrid[qi]^2
        # y = χ.*(wnGrid.*wnGrid)/qGrid[qi]^2
        # y0 = χ0 .* (wnGrid .* wnGrid) / qGrid[qi]^2
        # y = χ .* (wnGrid .* wnGrid) / qGrid[qi]^2

        NF = -Polarization.Polarization0_FiniteTemp(0.0, 0, param) * param.spin
        fxc = @. (1.0 / χ0 - 1 / χ) * NF / wnGrid / wnGrid * (para.qGrid[qi])^2
        # fxc = @. (1.0 / χ0 / 0.955 - 1 / χ) * NF

        filename = "Plot/chi_beta$(beta)_rs$(rs)_lam$(λ)_o$(order)_q$(qi-1).dat1"
        f = open(filename, "w")

        println("# q/kF: $(para.qGrid[qi]/kF)")
        @printf(f, "# %10s\t%10s\t%10s\t%10s\t%10s\n", "n", "ωn/EF", "Π", "Π0", "fxc")
        for (idn, n) in enumerate(nGrid)
            @printf(f, "%8i\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n", n, wnGrid[idn] / EF, χ[idn], χ0[idn], fxc[idn])
        end

        @printf("%10s\t%10s\t%10s\t%10s\t%10s\n", "n", "ωn/EF", "Π", "Π0", "fxc")
        for (idn, n) in enumerate(nGrid)
            @printf("%8i\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n", n, wnGrid[idn] / EF, χ[idn], χ0[idn], fxc[idn])
        end
        #plot(wnGrid, y0, curveconf = "w lp lw 1 lc '#08F7FE' pt 7 t 'χ_0 '")
        # plot!(f, wnGrid / EF, fxc, curveconf = "w lp lw 1 lc '#FFE64D' pt 7 t 'fxc'")
        #push!(p,plot(wnGrid, y0, curveconf = "w lp lw 1 lc '#FFE64D' pt 7 t 'χ_0'"))
    end
    #print(p)
    #plot(p)
    #save(term = "pdf",output= "myfigure.pdf",size = "1280,900",linewidth = 1)
    # display(f)
end

run()
