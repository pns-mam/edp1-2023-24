![PNS](http://caillau.perso.math.cnrs.fr/logo-pns.png)
## MAM4 - EDP1
# TP 2 - Équation de la chaleur

On souhaite simuler l'évolution de la température dans une pièce en forme de L
- dont le chauffage se situe dans la partie supérieure,
- dont les parois sont isolantes,
- à l'exception d'une fenêtre située sur le mur du bas.

On rappelle que cette température suit l'EDP de la chaleur 2D :

$$ \begin{aligned}
  & \frac{\partial u}{\partial t}(x,y,t) - \nu\Delta u(x,y,t) = f(x,y),\quad t > 0,\quad (x,y) \in \Omega,\\
  & u(x,y,0) = u_0(x,y),
\end{aligned} $$

avec des conditions aux limites appropriées sur le bord du domaine.

## Exercice 1

Compléter le code ci-dessous de façon à prendre en compte :
- le chauffage dans la zone $[0.6,0.9] \times [0.8,1]$,
- le fait que les parois sont isolantes (Neumann homogène),
- à l'exception de la fenêtre située entre $x = 0.25$ et $x = 0.75$ sur le mur sud (Dirichlet inhomogène).

NB. On pourra déplacer le chauffage pour tester différentes configurations.

```julia
using SparseArrays, Random, LinearAlgebra, Plots #; gr()

# Discrete 2D Laplacian - works with arbitrary ordering of interior nodes

function laplacian(G; Δx = 1.0, Δy = 1.0)
    ij = findall(G .≠ 0)
    K = length(ij)
    p = sortperm(G[ij])
    ij = ij[p] # ij is reordered so that G[ij[k]] = k
    Δ = spzeros(K, K)
    for k ∈ 1:K
        Δ[k, k] = -2 / Δx^2 - 2 / Δy^2
        ijₖ = ij[k]
        kW = G[ijₖ + CartesianIndex(0, -1)]
        if kW ≠ 0 Δ[k, kW] = 1 / Δx^2 end
        kE = G[ijₖ + CartesianIndex(0,  1)]
        if kE ≠ 0 Δ[k, kE] = 1 / Δx^2 end
        kS = G[ijₖ + CartesianIndex(-1, 0)]
        if kS ≠ 0 Δ[k, kS] = 1 / Δy^2 end
        kN = G[ijₖ + CartesianIndex( 1, 0)]
        if kN ≠ 0 Δ[k, kN] = 1 / Δy^2 end
    end
    return Δ
end

# Mesh generation
function meshgrid(x, y)
    n = length(x)
    m = length(y)
    xx = zeros(1, n)
    xx[:] = x
    yy = zeros(m, 1)
    yy[:] = y
    X = ones(m, 1) * xx
    Y = yy * ones(1, n)
    return (X, Y)
end

# Grid generation

ν = 1
Nx = 50
Ny = Nx
Nt = 2000
Nplot = 10
σ = 2 # CFL
x = range(0, 1, length = Nx + 1) # Ω is L-shaped included in (0,1)x(0,1)
y = range(0, 1, length = Ny + 1)
uw = -1 # window temperature
f = 250 # heating
Δx = x[2] - x[1]
Δy = y[2] - y[1]
Δt = σ / ν / (1 / Δx^2 + 1 / Δy^2)
tf = Nt  * Δt
(X, Y) = meshgrid(x, y)
ij = findall( (X .< 1) .& (Y .> 0) .& ( ( (X .> 0) .& (Y .< 0.5) ) .| ( (X .≥ 0.5) .& (Y .< 1) ) ) )
K = length(ij)
G = zeros(Int, size(X))
G[ij] = 1:K # columnwise ordering of interior nodes
spy(G[end:-1:1, :])

# Assembly

A = -ν * laplacian(G, Δx = Δx, Δy = Δy) # A ≃ -νΔ
b = zeros(K)
b[G[ (X .> .6) .& (X .< .9) .& (Y .> 0.3) .& (Y .< 0.7) ]] .= f

for k ∈ 1:K # Neumann everywhere...
    ijₖ = ij[k]
    kW = G[ijₖ + CartesianIndex(0, -1)]
    if kW == 0 A[k, k] = A[k, k] - 1 / Δx^2 end
    kE = G[ijₖ + CartesianIndex(0, 1)]
    if kE == 0 A[k, k] = A[k, k] - 1 / Δx^2 end
    kS = G[ijₖ + CartesianIndex(-1, 0)]
    if kS == 0
        if (X[ijₖ] ≥ 0.25) && (X[ijₖ] ≤ 0.75) # ... but on the window
            b[k] = b[k] + ν / Δy^2 * uw
        else
            A[k, k] = A[k, k] - 1 / Δy^2            
        end
    end
    kN = G[ijₖ + CartesianIndex(1, 0)]
    if kN == 0 A[k, k] = A[k, k] - 1 / Δy^2 end
end

# Solve for -νΔu = f (Poisson / stationary solution)

F = lu(A) # LU factorisation of sparse A (not assuming any symmetry)
u = F \ b
U = zeros(size(G))
U[ij] = u[G[ij]]
display( spy(A) )

# Plot solution

contour(x, y, U, fill=true, aspect_ratio=:equal)
```

## Exercice 2

Calculer l'évolution temporelle de la température à l'aide du schéma de Crank-Nicolson 2D :

$$ \begin{aligned}
  & u^{n+1}_k - (\sigma_x/2)(u^{n+1}_{kW}-2u^{n+1}_k+u^{n+1}_{kE})
              - (\sigma_y/2)(u^{n+1}_{kN}-2u^{n+1}_k+u^{n+1}_{kS})\\
\end{aligned} $$

$$ \begin{aligned}
  & u^{n+1}_k - (\sigma_x/2)(u^{n+1}_{kW}-2u^{n+1}_k+u^{n+1}_{kE})
              - (\sigma_y/2)(u^{n+1}_{kN}-2u^{n+1}_k+u^{n+1}_{kS})\\
  & = u^n_k   + (\sigma_x/2)(u^n_{kW}-2u^n_k+u^{n+1}_{kE})
              + (\sigma_y/2)(u^n_{kN}-2u^n_k+u^n_{kS}) + \Delta t \cdot f_k,
\end{aligned} $$

où $\sigma_x=\nu\Delta t/\Delta x^2$ et $\sigma_y=\nu\Delta t/\Delta y^2$.

```julia
# Evolution, Crank-Nicolson

u = uw * ones(K) # u0 = uniform initial temperature (consistent with Dirichlet) 
D = spdiagm(0 => 1.0 * ones(K)) + Δt/2 * A
E = spdiagm(0 => 1.0 * ones(K)) - Δt/2 * A
F = lu(D)

for i ∈ 1:Nt
    u = F \ (E * u + Δt * b)
    U[ij] = u[G[ij]]
    mod(i, Int(floor(Nt / Nplot))) == 0 && display( contour(x, y, U, fill=true, aspect_ratio=:equal) ) # contour plot
end
```