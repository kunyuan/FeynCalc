using Lehmann, LinearAlgebra, Printf, Parameters, DelimitedFiles
using ElectronGas

const N = 1000
const Euv, beta, rs = 100, 100, 1.0
const rtol = 1e-8
const EPS = 1e-10

if abspath(PROGRAM_FILE) == @__FILE__
    para = Parameter.rydbergUnit(1 / beta, rs)
    EF = para.EF
    dlr = DLRGrid(EF * Euv, beta / EF, rtol, false, :ph)

    # grid = readdlm("ph_10000_1e-8.dlr", comments = true, comment_char = '#')
    # τgrid = grid[:, 3]

    τgrid = dlr.τ / para.β
    # τgrid = [0.0, 0.5]
    len = length(τgrid)
    grid = zeros(Float64, (len - 1) * N + 1)
    x = range(1, (len - 1) * N + 1)

    local ti = 2
    grid[1] = τgrid[1]
    for (τi, τ) in enumerate(τgrid)
        if τi < len
            δτ = τgrid[τi+1] - τ
            for j in 1:N
                grid[ti] = τ + j * δτ / N
                ti = ti + 1
            end
        end
    end
    # grid[1] = EPS
    # grid[end, :] = [ti τgrid[end]]
    open("ph_tau.dlr", "w") do io
        writedlm(io, zip(x, grid))
    end
end