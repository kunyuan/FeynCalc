# Calculate the susceptibility vs. matsubara frequency 
# using QuantumStatistics, LinearAlgebra, Printf, Parameters, DelimitedFiles
using Lehmann, LinearAlgebra, Printf, Parameters, DelimitedFiles
using Roots, Polylogarithms, FastGaussQuadrature
using Gaston, LaTeXStrings, ElectronGas

const Qsize = 6

function chemical_potential(beta)
    f(β, μ) = real(polylog(3 / 2, -exp(β * μ))) + 4 / 3 / π^0.5 * (β)^(3 / 2)
    g(μ) = f(beta, μ)
    return find_zero(g, (-10000, 1))
end

@with_kw struct Para
    beta::Float64
    rs::Float64
    λ::Float64
    order::Int

    kF::Float64
    EF::Float64 = kF^2
    NF::Float64 = kF / 2 / π^2
    β::Float64 = beta / EF
    qGrid::Vector{Float64}
    nsize::Int = 20
    nGrid::Vector{Int64} = [n for n in 0:nsize-1]
    μ = chemical_potential(beta) * EF
    # χdlr = DLRGrid(EF * 10, beta / EF, 1e-3, false, :ph)
end

function readPara(line)
    # paras = readdlm(fname, '\n')
    # param = parse.(Float64, split.(line[1, :])[1])
    param = parse.(Float64, split.(line)[:])
    order, beta, rs, λ = param[1], param[2], param[3], param[5]
    minQ, maxQ = param[6], param[7]
    kF = (9π / 4)^(1 / 3) / rs
    qGrid = [q for q in LinRange(minQ, maxQ - (maxQ - minQ) / Qsize, Qsize)]

    return Para(beta=beta, rs=rs, λ=λ, order=order, kF=kF, qGrid=qGrid)
end

function susceptibility(para, QBin)
    @unpack beta, β, rs, λ, order, kF, EF, qGrid, nsize, nGrid, μ = para

    χ0 = zeros(Float64, nsize)
    # χ = similar(χ0)
    q = qGrid[QBin] * kF

    filename = "Data_Polar/polar_beta$(beta)_rs$(rs)_lam$(λ)_o$(order)_q$(QBin-1)_Uni.dat2"
    print("Loading ", filename, "\n")
    polar_tau = readdlm(filename)
    χ_tau::Vector{Float64} = polar_tau[6:end, 2]
    τGrid::Vector{Float64} = polar_tau[6:end, 1]

    rtol = polar_tau[6, 3] / polar_tau[6, 2]
    χdlr = DLRGrid(EF * 10, beta / EF, rtol, false, :ph)

    @printf("ExtMom: %f, DLR nsize: %d\n", q, χdlr.size)
    χ = tau2matfreq(χdlr, χ_tau, para.nGrid, τGrid)
    #χ = DLR.tau2matfreq(:corr, χ_tau, para.χdlr, para.nGrid, axis=1)
    #print(tau2matfreq(para.χdlr,χ_tau))
    print(" n χ_0 χ \n")
    param = Parameter.rydbergUnit(1 / beta, rs, μ=μ)
    for (idn, n) in enumerate(nGrid)
        χ0[idn] = Polarization.Polarization0_FiniteTemp(q, n, param) * 2
        @printf("%d %10.8f %10.8f \n", n, -χ0[idn], χ[idn])
    end
    return -χ0, real(χ)
end

function run(line)
    para = readPara(line)
    @unpack beta, rs, λ, order, qGrid, nGrid, kF, β, EF, NF = para

    p = Gaston.Figure[]
    Fig = Gaston.Figure[]
    wnGrid = nGrid * 2π / β
    for qi in 2:Qsize
        println("ω_0 = $(2π/β),  vF*q = $(qGrid[qi]*EF*2)")
        χ0, χ = susceptibility(para, qi)
        # y0 = χ0 .* (wnGrid .* wnGrid) / (qGrid[qi] * kF)^2
        # y = χ .* (wnGrid .* wnGrid) / (qGrid[qi] * kF)^2
        # filename = "Plot/pi_beta$(beta)_rs$(rs)_lam$(λ)_o$(order)_q$(qi-1).dat"
        # f = open(filename, "w")
        # f_xc = zeros(Float64, para.nsize)
        f_xc = 1 ./ χ - 1 ./ χ0
        # @printf(f, "%d\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n", n, wnGrid[idn], y[idn], y0[idn], χ[idn], χ0[idn], f_xc[idn])

        x = wnGrid / EF
        y = f_xc * NF
        # y = f_xc * NF * (qGrid[qi] * kF)^2 ./ (wnGrid .* wnGrid)
        # Figmid =
        if qi == 2
            plot(x, y, curveconf="tit 'q=$(qGrid[qi]) kF' w lp lw 1 lc '#08F7FE'",
                Axes(# xrange = (0, 0.8),
                    # yrange = (-4, 1.2),
                    yrange=(-0.2, -0.05),
                    ylabel="'f_{xc} N_F'",
                    # ylabel = "'f_{xc} N_F (q/ω_n)^2'",
                    xlabel="'ω_n/E_F'",
                    key="b l",
                    title="'beta=$beta, r_s=$rs, λ=$λ'")
            )
        else
            plot!(x, y, curveconf="tit 'q=$(qGrid[qi]) kF' w lp lw 1")
        end
        # push!(Fig, Figmid)
        # save(term = "pdf", output = "./Plot/fxc_beta$(beta)_rs$(rs)_q$(qi)_qm.pdf", linewidth = 1)
    end
    # plot([Fig[1] Fig[2]; Fig[3] Fig[4]])
    save(term="pdf", output="./Plot/fxc_beta$(beta)_rs$(rs)_lam$(λ)_o$(order)_qm_Uni.pdf", linewidth=1)
    # save(term = "pdf", output = "./Plot/fxcRe_beta$(beta)_rs$(rs)_lam$λ.pdf", linewidth = 1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    for line in eachline("inlist")
        length(line) == 0 && break
        run(line)
    end
end