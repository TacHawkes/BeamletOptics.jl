struct StepIndexFiber{T}
    radius::T
    n_core::T
    n_clad::T
end

"""
    mode_field_radius(fiber, λ)

Estimate the fundamental mode field radius of `fiber` at wavelength `λ`
using the empirical formula by Petermann.
"""
function mode_field_radius(f::StepIndexFiber{T}, λ::Real) where T
    V = 2π/λ * f.radius * sqrt(f.n_core^2 - f.n_clad^2)
    return f.radius * (0.65 + 1.619/V^(3/2) + 2.879/V^6)
end

mode_field_diameter(f::StepIndexFiber, λ) = 2*mode_field_radius(f, λ)

"""
    fundamental_mode(fiber, λ, r)

Return the normalized fundamental mode field of `fiber` evaluated at the
radial coordinate(s) `r`.
"""
function fundamental_mode(f::StepIndexFiber, λ, r)
    w = mode_field_radius(f, λ)
    A = sqrt(2/(π*w^2))
    return A .* exp.(-(r.^2)./(w^2))
end

"""
    electric_field(pd::Photodetector)

Return the complex electric field accumulated on `pd`.
"""
electric_field(pd::Photodetector) = pd.field

function electric_field(psf::PSFDetector{T};
        n::Int=100,
        crop_factor::Real=1,
        center::Symbol=:centroid,
        x_min = Inf,
        x_max = Inf,
        z_min = Inf,
        z_max = Inf,
        x0_shift::Real=0,
        z0_shift::Real=0) where T
    _x_min, _x_max, _z_min, _z_max = calc_local_lims(psf; crop_factor=crop_factor, center=center)
    if x_min != Inf && x_max != Inf
        _x_min = x_min; _x_max = x_max
    end
    if z_min != Inf && z_max != Inf
        _z_min = z_min; _z_max = z_max
    end
    xs = LinRange(_x_min, _x_max, n) .+ x0_shift
    zs = LinRange(_z_min, _z_max, n) .+ z0_shift
    field = zeros(Complex{T}, n, n)
    orient = orientation(psf)
    @views e1, e2 = Point3(orient[:, 1]), Point3(orient[:, 3])
    origin_pd = position(psf)
    Threads.@threads for j in eachindex(zs)
        z = zs[j]
        @inbounds for i in eachindex(xs)
            x = xs[i]
            p = origin_pd + x * e1 + z * e2
            acc = zero(complex(T))
            @inbounds @simd for h in psf.data
                l = dot(p - position(h), direction(h))
                acc += projection_factor(h) * cis(wavenumber(h) * (optical_path_length(h) + l))
            end
            field[i,j] = acc
        end
    end
    return xs, zs, field
end

"""
    fiber_coupling(pd::Photodetector, fiber, λ)

Compute the coupling efficiency between the field measured by `pd` and the
fundamental mode of `fiber` at wavelength `λ`.
"""
function fiber_coupling(pd::Photodetector, f::StepIndexFiber, λ)
    field = electric_field(pd)
    xs = pd.x; ys = pd.y
    dx = step(xs); dy = step(ys)
    r = [sqrt(x^2 + y^2) for x in xs, y in ys]
    mode = fundamental_mode(f, λ, r)
    dA = dx*dy
    num = abs(sum(conj.(field) .* mode) * dA)^2
    denom = (sum(abs2, field) * dA) * (sum(abs2, mode) * dA)
    return num / denom
end

function fiber_coupling(psf::PSFDetector, f::StepIndexFiber, λ; kwargs...)
    xs, zs, field = electric_field(psf; kwargs...)
    dx = step(xs); dz = step(zs)
    r = [sqrt(x^2 + z^2) for x in xs, z in zs]
    mode = fundamental_mode(f, λ, r)
    dA = dx*dz
    num = abs(sum(conj.(field) .* mode) * dA)^2
    denom = (sum(abs2, field) * dA) * (sum(abs2, mode) * dA)
    return num / denom
end

const coupling_efficiency = fiber_coupling

