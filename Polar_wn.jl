# Calculate the susceptibility vs. matsubara frequency 
# using QuantumStatistics, LinearAlgebra, Printf, Parameters, DelimitedFiles
using Lehmann, LinearAlgebra, Printf, Parameters, DelimitedFiles
using Roots, Polylogarithms, FastGaussQuadrature
using Gaston, LaTeXStrings, ElectronGas
#using QuantumStatistics

# include("parameter.jl")
const Qsize = 5

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
    β::Float64 = beta / EF
    # qGrid::Vector{Float64} = [q for q in LinRange(0.0, 0.08 * kF, Qsize)]
    qGrid::Vector{Float64}
    nsize::Int = 11
    nGrid::Vector{Int64} = [n for n in 0:nsize-1]
    μ = chemical_potential(beta) * EF
    χdlr = DLRGrid(EF * 100, beta / EF, 1e-8, false, :ph)
end

function readPara(fname)
    line = readdlm(fname, '\n')
    param = parse.(Float64, split.(line[1, :])[1])
    order, beta, rs, λ = param[1], param[2], param[3], param[5]
    minQ, maxQ = param[5], param[6]
    kF = (9π / 4)^(1 / 3) / rs
    qGrid = [q for q in LinRange(0.0, (maxQ - (maxQ - minQ) / Qsize) * kF, Qsize)]

    return Para(beta = beta, rs = rs, λ = λ, order = order, kF = kF, qGrid = qGrid)
end

function susceptibility(para, QBin)
    @unpack beta, β, rs, λ, order, qGrid, nsize, nGrid, μ, χdlr = para

    χ0 = zeros(Float64, nsize)
    # χ = similar(χ0)
    q = qGrid[QBin]
    @printf("ExtMom: %f, DLR nsize: %d\n", q, χdlr.size)

    filename = "Data_Polar/polar_beta$(beta)_rs$(rs)_lam$(λ)_o$(order)_q$(QBin-1).dat1"
    print("Loading ", filename, "\n")
    polar_tau = readdlm(filename)
    χ_tau::Vector{Float64} = polar_tau[6:end, 2]
    #print(χ_tau)    
    #print(χ)
    χ = tau2matfreq(χdlr, χ_tau, para.nGrid)
    #χ = DLR.tau2matfreq(:corr, χ_tau, para.χdlr, para.nGrid, axis=1)
    #print(tau2matfreq(para.χdlr,χ_tau))
    print(" n χ_0 χ \n")
    param = Parameter.rydbergUnit(1 / β, rs, μ = μ)
    for (idn, n) in enumerate(nGrid)
        χ0[idn] = Polarization.Polarization0_FiniteTemp(q, n, param) * 2
        @printf("%d %10.8f %10.8f \n", n, -χ0[idn], χ[idn])
    end
    return -χ0, real(χ), param
end

function run()
    para = readPara("inlist")

    @unpack beta, rs, qGrid, nGrid, kF, β = para
    p = Gaston.Figure[]
    Fig = Gaston.Figure[]
    wnGrid = nGrid * 2π / β
    for qi in 2:Qsize
        χ0, χ = susceptibility(para, qi)
        y0 = χ0 .* (wnGrid .* wnGrid) / qGrid[qi]^2
        y = χ .* (wnGrid .* wnGrid) / qGrid[qi]^2
        # filename = "Plot/pi_beta$(beta)_rs$(rs)_lam$(λ)_o$(order)_q$(qi-1).dat"
        # f = open(filename, "w")
        f_xc = zeros(Float64, para.nsize)
        for (idn, n) in enumerate(nGrid)
            f_xc[idn] = 1 ./ χ[idn] - 1 ./ χ0[idn]
            # @printf(f, "%d\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n", n, wnGrid[idn], y[idn], y0[idn], χ[idn], χ0[idn], f_xc[idn])
        end
        #println(f_xc)
        plot(wnGrid, y0, curveconf = "w lp lw 1 lc '#08F7FE' pt 7 t 'χ_0'",
            Axes(object = "rectangle from screen 0,0 to screen 1,1 behind fc 'black' fs solid noborder",
                border = "lw 1 lc 'white'",
                xtics = "textcolor rgb 'white'",
                ytics = "textcolor rgb 'white'",
                title = "'beta=$(beta), rs=$(rs), q=$(qGrid[qi]/kF) kF' textcolor 'white'",
                ylabel = "'χω_n^2/q^2' textcolor 'white'",
                xlabel = "'ω_n' textcolor 'white'",
                grid = "ls 1 lc '#2A3459'",
                key = "t r textcolor 'white'"),
            handle = qi
        )
        pmid = plot!(wnGrid, y, curveconf = "w lp lw 1 lc '#FFE64D' pt 7 t 'Π'", handle = qi)
        push!(p, pmid)
        save(term = "pdf", output = "./Plot/polar_wn_beta$(beta)_rs$(rs)_q$qi.pdf", linewidth = 1)

        Figmid = plot(wnGrid, f_xc, curveconf = "w lp lw 1 lc '#08F7FE' pt 7 ",
            Axes(object = "rectangle from screen 0,0 to screen 1,1 behind fc 'black' fs solid noborder",
                border = "lw 1 lc 'white'",
                xtics = "textcolor rgb 'white'",
                ytics = "textcolor rgb 'white'",
                title = "'beta=$(beta), rs=$(rs), q=$(qGrid[qi]/kF) kF' textcolor 'white'",
                ylabel = "'v_qG_-' textcolor 'white'",
                xlabel = "'ω_n' textcolor 'white'",
                grid = "ls 1 lc '#2A3459'",
                key = "t r textcolor 'white'"),
            handle = qi + 5
        )
        push!(Fig, Figmid)
        save(term = "pdf", output = "./Plot/fxc_beta$(beta)_rs$(rs)_q$qi.pdf", linewidth = 1)
    end
    # plot([p[1] p[2]; p[3] p[4]], handle = 11)
    plot([p[1] p[2]; p[3] p[4]])
    save(term = "pdf", output = "./Plot/polar_wn_beta$(beta)_rs$(rs).pdf")
    # plot([Fig[1] Fig[2]; Fig[3] Fig[4]], handle = 12)
    plot([Fig[1] Fig[2]; Fig[3] Fig[4]])
    save(term = "pdf", output = "./Plot/fxc_beta$(beta)_rs$(rs).pdf", linewidth = 1)

    # plot([Fig[1]])
    # save(term = "pdf", output = "./Plot/fxc_beta$(beta)_rs$(rs)_q2.pdf", linewidth = 1)

end

run()