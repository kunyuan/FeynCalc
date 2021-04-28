using Cuba, QuantumStatistics, Trapz, Roots

###### constants ###########
const e0 = sqrt(2)  # electric charge
const me = 0.5  # electron mass
const dim = 3    # dimension (D=2 or 3, doesn't work for other D!!!)
const spin = 2  # number of spins

const rs = 1.0  
const kF = (dim == 3) ? (9π / (2spin))^(1 / 3) / rs : sqrt(4 / spin) / rs
const EF = kF^2 / (2me)
const β = 1.0 / kF^2
const mass2 = 1.0
const maxK = 6kF

println("rs=$rs, β=$β, kF=$kF, EF=$EF, mass2=$mass2")

const kgrid = Grid.fermiK(kF, maxK, 0.1kF, 512)  # external K grid for sigma

@inline function linear1D(data, xgrid, x)
    xarray = xgrid.grid
    if (x <= xarray[firstindex(xgrid)]) # the first k bin is only slightly larger than 0. Therefore, any k smaller than the frist bin can be regarded as the first bin
        return data[1]
    end
    if (x >= xarray[lastindex(xgrid)])
        return 0.0
    end
    xi0 = floor(xgrid, x)
    xi1 = xi0 + 1
    dx0, dx1 = x - xarray[xi0], xarray[xi1] - x

    d0, d1 = data[xi0], data[xi1]

    return (data[xi0] * dx1 + data[xi1] * dx0) / (dx0 + dx1)
end

@inline function fockT0()
    fock = similar(kgrid.grid)
    l = sqrt(mass2)
    for (ki, k) in enumerate(kgrid.grid)
        if ki == 1
            k = 1.0e-3
        end
        fock[ki] = 1.0 + l / kF * atan((k - kF) / l);
        fock[ki] -= l / kF * atan((k + kF) / l);
        fock[ki] -= (l * l - k * k + kF * kF) / 4.0 / k / kF *
                log((l * l + (k - kF) * (k - kF)) / (l * l + (k + kF) * (k + kF)));
        fock[ki] *= (-2.0 * kF) / π;
    end

    k = kF
    dμ = 1.0 + l / kF * atan((k - kF) / l);
    dμ -= l / kF * atan((k + kF) / l);
    dμ -= (l * l - k * k + kF * kF) / 4.0 / k / kF *
            log((l * l + (k - kF) * (k - kF)) / (l * l + (k + kF) * (k + kF)));
    dμ *= (2.0 * kF) / π;

return fock, dμ
end

function FockStatic(dim::Int, k::T, β::T, me::T, mass2::T, spin, sigmaK, dμ::T) where {T <: AbstractFloat}
    if k < 0.0
        k = -k
    end

    if k / kF < 1.0e-6
        k = 1.0e-6 * kF
    end

    function fock(q)
        phase = T(1.0)
            if dim == 3
            phase *= 1 / (π * k)
        else
            error("not implemented")
        end
        sigma = linear1D(sigmaK, kgrid, q)
        ϵ = β * (q^2 - kF^2 + sigma + dμ) / (2me)

        p = - phase * Spectral.fermiDirac(ϵ) * q * log(((k + q)^2 + mass2^2) / ((k - q)^2 + mass2^2))

            if isnan(p)
            println("warning: integrand at ω=$ω, q=$q, k=$k is NaN!")
        end
        # println(p)
        return p
    end

    function integrand(x, f)
        # x[1]:k
        f[1] = fock(x[1] / (1 - x[1])) / (1 - x[1])^2
    end

    result, err = Cuba.cuhre(integrand, 2, 1, rtol=1.0e-6)
    # result, err = Cuba.vegas(integrand, 1, 1, rtol=rtol)
    return result[1], err[1]
end

function density(sigmaK, dμ)
    nk = similar(sigmaK)
    for (ki, k) in enumerate(kgrid.grid)
        ϵ = k^2 - kF^2 + sigmaK[ki] + dμ
        # ϵ = k^2 - kF^2 + dμ
        nk[ki] = Spectral.fermiDirac(β * ϵ) * 4π * k^2 / (2π)^3 * spin
    end
    return trapz(collect(kgrid.grid), nk)
end
    
    function fockK(dμ)
    oldsigmaK, _ = fockT0()
    newsigmaK = similar(kgrid.grid)
    for i in 1:10
        for (ki, k) in enumerate(kgrid.grid)
            newsigmaK[ki] = FockStatic(dim, k, β, me, mass2, spin, oldsigmaK, dμ)[1]
        end
        error = maximum(abs.(newsigmaK - oldsigmaK))
        println("iteration $i -> max error: ", error)
        if error < 1.0e-4
            break
        end
        oldsigmaK = copy(newsigmaK)
    end
    return newsigmaK
end
    
function diff(dμ)
    sigma = fockK(dμ)
    n = density(sigma, dμ)
    n0 = 1 / (4 / 3 * π * rs^3) # exact density
    println("density: $n, expected $n0")
    return n - n0
end

if abspath(PROGRAM_FILE) == @__FILE__
    # using Gaston
    # fock = fockK(EF)
    # for (ki, k) in enumerate(kgrid.grid)
    #     println("$(k / kF)   $(fock[ki])")
    # end

    sigmaK, dμ0 = fockT0()
    # sigmaK = fockK(dμ0)
    # n0 = 1 / (4 / 3 * π * rs^3) # exact density

    # println(density(sigmaK, dμ0), " expected density = $n0")

    dμ = find_zero(diff, (0.0, maxK), Bisection(), rtol=1.0e-4)
    println("finite T : $dμ, zero T: $dμ0")

    sigmaK = fockK(dμ)

    for (ki, k) in enumerate(kgrid.grid)
        println("$(k / kF)   $(sigmaK[ki])  $dμ")
    end

    # for dμ in LinRange(dμ0 - kF / 3, dμ0 + kF / 3, 10)
    #     println(diff(dμ))
    # end

    # for (ki, k) in enumerate(kgrid.grid)
    #     println("$(k / kF)   $(sigmaK[ki])")
    # end
end