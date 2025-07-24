# Fiber coupling example

This example demonstrates how to estimate the coupling efficiency of a Gaussian
beam into a single--mode fiber.  The electric field at the fiber facet is
sampled with a [`Photodetector`](@ref) and compared to the fundamental mode of a
[`StepIndexFiber`](@ref).

```@example fiber_coupling
using BeamletOptics

λ = 1064e-9
fiber = StepIndexFiber(4.1e-6, 1.46, 1.455)

w = mode_field_radius(fiber, λ)
beam = GaussianBeamlet([0, -0.01, 0], [0, 1, 0], λ, w)

pd = Photodetector(25e-6, 200)
system = System([pd])
solve_system!(system, beam)

η = fiber_coupling(pd, fiber, λ)
println("coupling = ", η)
```
