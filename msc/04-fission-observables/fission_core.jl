# Quantities shared by the scripts of this directory: the compound system, the
# isobaric charge distribution, Q values and separation energies from the mass
# table, and the reduction of the yield matrix to the (A_H, TKE) cells that were
# measured.

include(joinpath(@__DIR__, "fission_data.jl"))

const DATA = joinpath(@__DIR__, "data")
const FIGURES = joinpath(@__DIR__, "figures")

"Mass number of the compound system ²³⁶U."
const A₀ = 236
"Proton number of the compound system."
const Z₀ = 92
"Charge polarisation ΔZ of the heavy fragment; the light fragment takes +0.5."
const ΔZ_POL = -0.5
"Rms width of the Gaussian isobaric charge distribution."
const σ_Z = 0.6
"Heavy-fragment mass range of the exercises."
const A_H_RANGE = 118:160

"Most probable heavy-fragment charge, Z_p(A_H) = Z_UCD(A_H) + ΔZ."
most_probable_charge(A_H) = Z₀ * A_H / A₀ + ΔZ_POL

"Unnormalised Gaussian weight of the heavy-fragment charge `Z_H` about Z_p(A_H)."
charge_weight(Z_H, A_H) = exp(-(Z_H - most_probable_charge(A_H))^2 / (2σ_Z^2))

"The three integer charges nearest Z_p(A_H)."
function nearest_three_charges(A_H)
    centre = round(Int, most_probable_charge(A_H))
    return (centre - 1):(centre + 1)
end

"Every charge for which the mass table lists the isobar `A`."
tabulated_charges(masses, A) = sort!([Z for (Z, a) in keys(masses) if a == A])

"""
    q_value(masses, Z_H, A_H) -> (Q, σQ) or nothing

Q value in MeV of the split ²³⁶U → (Z_H, A_H) + (Z₀ − Z_H, A₀ − A_H) from the
mass excesses, with the three mass uncertainties in quadrature. `nothing` when
either fragment is absent from the table.
"""
function q_value(masses, Z_H, A_H)
    heavy = get(masses, (Z_H, A_H), nothing)
    light = get(masses, (Z₀ - Z_H, A₀ - A_H), nothing)
    (heavy === nothing || light === nothing) && return nothing
    compound = masses[(Z₀, A₀)]
    Q = (compound[1] - heavy[1] - light[1]) / 1000
    return Q, sqrt(compound[2]^2 + heavy[2]^2 + light[2]^2) / 1000
end

"""
    separation_energy(masses, Z, A) -> (S_n, σ) or nothing

Neutron separation energy of the nuclide `(Z, A)` in MeV.
"""
function separation_energy(masses, Z, A)
    parent = get(masses, (Z, A), nothing)
    daughter = get(masses, (Z, A - 1), nothing)
    (parent === nothing || daughter === nothing) && return nothing
    neutron = masses[(0, 1)]
    S = (daughter[1] + neutron[1] - parent[1]) / 1000
    return S, sqrt(parent[2]^2 + daughter[2]^2 + neutron[2]^2) / 1000
end

"""
    charge_average(f, A_H, charges) -> (mean, σ) or nothing

Average of `f(Z_H) -> (value, σ)` over `charges` with the Gaussian weights of the
isobaric charge distribution; charges for which `f` returns `nothing` are left
out. The uncertainty combines the `σ` of the terms in quadrature, as the
originals did.
"""
function charge_average(f, A_H, charges)
    num = 0.0
    den = 0.0
    var = 0.0
    for Z_H in charges
        term = f(Z_H)
        term === nothing && continue
        w = charge_weight(Z_H, A_H)
        num += w * term[1]
        var += (w * term[2])^2
        den += w
    end
    return den > 0 ? (num / den, sqrt(var) / den) : nothing
end

"One measured cell of the yield matrix: Y(A_H, TKE) and its error, in %."
struct YieldCell
    A_H::Int
    TKE::Int
    Y::Float64
    σY::Float64
end

"""
    yield_cells(rows) -> Vector{YieldCell}

The (A_H, TKE) cells of the yield matrix. The five charge rows of a cell are one
measured yield split by fixed fractions — they share one relative error — so
within a cell the errors add linearly.
"""
function yield_cells(rows)
    acc = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    for r in rows
        Y, σ = get(acc, (r.A_H, r.TKE), (0.0, 0.0))
        acc[(r.A_H, r.TKE)] = (Y + r.Y, σ + r.σY)
    end
    return [YieldCell(k[1], k[2], acc[k]...) for k in sort!(collect(keys(acc)))]
end

"""
    marginal(rows, key) -> (keys, Y, σY)

Yield summed over every row sharing `key(row)`. Rows of one (A_H, TKE) cell are
fully correlated and add linearly; distinct cells are independent and add in
quadrature.
"""
function marginal(rows, key)
    Y = Dict{Int,Float64}()
    σ_cell = Dict{Tuple{Int,Int,Int},Float64}()
    for r in rows
        k = key(r)
        Y[k] = get(Y, k, 0.0) + r.Y
        σ_cell[(k, r.A_H, r.TKE)] = get(σ_cell, (k, r.A_H, r.TKE), 0.0) + r.σY
    end
    var = Dict{Int,Float64}()
    for ((k, _, _), σ) in σ_cell
        var[k] = get(var, k, 0.0) + σ^2
    end
    ks = sort!(collect(keys(Y)))
    return ks, [Y[k] for k in ks], [sqrt(var[k]) for k in ks]
end

"""
    yield_weighted_mean(q, Y, σY) -> (mean, σ)

Mean of `q` weighted by the independent yields `Y`, and the uncertainty the
yield errors give it,

    ⟨q⟩ = Σ qᵢYᵢ / ΣY,        σ² = Σ [(qᵢ − ⟨q⟩) σYᵢ / ΣY]²,

which is the first-order propagation of ∂⟨q⟩/∂Yᵢ = (qᵢ − ⟨q⟩)/ΣY.
"""
function yield_weighted_mean(q, Y, σY)
    ΣY = sum(Y)
    mean = sum(q .* Y) / ΣY
    return mean, sqrt(sum(abs2, (q .- mean) .* σY)) / ΣY
end

"⟨TKE⟩ [MeV] and yield [%] of every heavy-fragment mass with a measured yield."
function tke_by_mass(cells)
    Y = Dict{Int,Float64}()
    T = Dict{Int,Float64}()
    for c in cells
        Y[c.A_H] = get(Y, c.A_H, 0.0) + c.Y
        T[c.A_H] = get(T, c.A_H, 0.0) + c.Y * c.TKE
    end
    return Dict(A => (TKE = T[A] / Y[A], Y = Y[A]) for A in keys(Y) if Y[A] > 0)
end

"""
    assert_close(name, value, reference; atol = 0, rtol = 0)

Throw unless `value ≈ reference` within the given tolerance.
"""
function assert_close(name, value, reference; atol = 0.0, rtol = 0.0)
    isapprox(value, reference; atol = atol, rtol = rtol) ||
        throw(ErrorException("$name: $value against $reference " *
                             "(atol = $atol, rtol = $rtol)"))
    return nothing
end
