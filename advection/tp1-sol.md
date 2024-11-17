![PNS](http://caillau.perso.math.cnrs.fr/logo-pns.png)
## MAM4 - EDP1
# TP 2 - Équation d'advection

## Exercice 1. Schéma explicite

On résout numériquement l'équation d'advection sur $\Omega := ]0,1[$ avec conditions aux limites périodiques ($V > 0$) :

$$ \left\lbrace \begin{align}
& \frac{\partial u}{\partial t}(t,x) + V\frac{\partial u}{\partial x}(t,x) = 0,\quad x \in \Omega,\quad t > 0,\\
& u(t,x) = u(t,x+1),\quad t \ge 0,\\
& u(0,x) = u_0(x),\quad x \in \Omega.
\end{align} \right. $$

On suppose la donnée initiale $u_0$ $1$-périodique. On vérifie aisément que la
solution exacte est donnée par $u(t,x) = u_0(x-Vt)$.

On cherche à approcher numériquement la solution par le schéma *décentré amont* (*cf.* $V>0$) suivant :

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + V\frac{u_j^n-u_{j-1}^n}{\Delta x} = 0,\quad J \geq 1,\quad n \geq 0, $$ 

où $u_j^n \simeq u(t_n,x_j)$, $t_n = n\Delta t$ et $x_j = j\Delta x$.
En réarrangeant les termes, on obtient

$$ u_j^{n+1} = u_j^n - \sigma(u_j^n-u_{j-1}^n) $$

où $\sigma = V\Delta t/\Delta x$ est connu sous le nom de
[*nombre de Courant*](https://fr.wikipedia.org/wiki/Nombre_de_Courant).
Pour approcher la condition initiale et la condition limite on écrit

$$ u_j^0 = u_0(x_j),\quad 0 \le j \le J,\quad n \ge 0. $$

Augmenter progressivement le paramètre $\sigma$ et observer le résultat. Quelle est la valeur critique ? Constater aussi en augmentant progressivement $N$ que la solution numérique est amortie au fil des itérations en temps (phénomène de diffusion numérique).

```julia
using LinearAlgebra, SparseArrays, Plots
  
# Parameters
V = 0.1                                            # advection speed
J = 1600                                           # space grid size
x = range(0, 1, length=J+1)                        # space grid
Δx = x[2] - x[1]                                   # space stepsize
σ = 1.1                                            # Courant number
Δt = σ * Δx / V                                    # time stepsize
N = 2400                                           # time grid size     
tf = N * Δt                                        # final time
function u0(x)
    y = mod(x, 1)
    return Float64((y > 0.2) & (y < 0.3))
end

# Explicit scheme
u = u0.(x)

for n ∈ 1:N
    u[2:end] = u[2:end] - σ * (u[2:end] - u[1:end-1])
    u[1] = u[end]
end

uexact = u0.(x .- V * tf)
err = u - uexact
u_plot = plot(x, u, xlabel="x", ylabel="u", color=:red, label="Finite differences", lw=6)
plot!(u_plot, x, uexact, xlabel="x", ylabel="u", color=:black, label="Solution", lw=2)
err_plot = plot(x, err, xlabel="x", ylabel="Error", legend=false)
display(plot(u_plot, err_plot, layout=(2, 1), size=(700, 700)))
println("Δx: ", Δx, "\t Δt:", Δt, "\t max error: ", maximum(abs.(err)))
```

## Exercice 2. Schéma de Lax-Wendroff
En suivant le modèle précédent, implémenter le schéma suivant :

$$ u_j^{n+1} = u_j^n - \frac{\sigma}{2}(u_{j+1}^n - u_{j-1}^n) + \frac{\sigma^2}{2}(u_{j+1}^n - 2u_j^n + u_{j-1}^n). $$

Dans ce nouveau script au tout début, à l'intérieur de la boucle en temps on va fixer la condition à un des bords (condition entrante) :

```julia
u[1] = u0(x[1] - V * n * Δt)
u[end] = u[1]
```

Repartir de $\sigma = 0.8$ et tester plusieurs possibilités du couple conditions initiale / schéma numérique. 
- Quelles conditions tirez-vous ?
- Augmenter progressivement $\sigma$ pour le schéma de Lax-Wendroff et observer. 
- Que peut-on dire de la diffusion numérique observée précédemment dans le cas du schéma décentré ?
- Le schéma de Lax-Wendroff est-il diffusif ?

```julia
# Explicit Lax-Wendroff scheme
u = u0.(x)

for n ∈ 1:N
    u[2:end-1] = u[2:end-1] - .5σ * (u[3:end] - u[1:end-2]) + .5σ^2 * (u[3:end] - 2u[2:end-1] + u[1:end-2])
    u[1] = u0(x[1] - V * n * Δt)
    u[end] = u[1]
end

uexact = u0.(x .- V * tf)
err = u - uexact
u_plot = plot(x, u, xlabel="x", ylabel="u", color=:red, label="Finite differences", lw=6)
plot!(u_plot, x, uexact, xlabel="x", ylabel="u", color=:black, label="Solution", lw=2)
err_plot = plot(x, err, xlabel="x", ylabel="Error", legend=false)
display(plot(u_plot, err_plot, layout=(2, 1), size=(700, 700)))
println("Δx: ", Δx, "\t Δt:", Δt, "\t max error: ", maximum(abs.(err)))
```

## Exercice 3. Schéma de Lax-Wendroff implicite
On va changer le schéma comme suit :

$$ u_j^{n+1} = u_j^n - \frac{\sigma}{2}(u_{j+1}^n - u_{j-1}^n)
             + \frac{\sigma^2}{2}(u_{j+1}^{n+1} - 2u_j^{n+1} + u_{j-1}^{n+1}). $$ 

En introduisant $W^n := (u_j^n - (\sigma/2)(u_{j+1}^n - u_{j-1}^n))_j$ et $U^n := (u_j^n)_j$,
le schéma s'écrit $AU^{n+1} = W^n$. On pourra utiliser la
condition de périodicité directement dans le schéma numérique afin d'éliminer $u_j^N$. On constate que le schéma est implicite : l'évaluation de $U^{n+1}$ à partir
de $U^n$ nécessite la résolution d'un système linéaire, et le coût
de l'itération sera plus élevé que dans le cas d'un
schéma explicite. Cependant le schéma n'est pas limité par la
valeur du nombre de Courant $\sigma$ ce qui permet d'utiliser des
pas de temps plus grands et donc de diminuer le nombre d'itérations associé.

Réaliser des expériences numériques en
faisant varier le nombre de Courant. Discuter les performances relatives des schémas implicite et explicite en termes de temps de calcul.

```julia
# Implicit Lax-Wendroff scheme
u = u0.(x)

A = spdiagm(-1 => -σ^2/2*ones(J), 0 => (1+σ^2)*ones(J+1), 1 => -σ^2/2*ones(J))
A[1, end] = -.5σ^2
A[end, 1] = -.5σ^2

F = cholesky(A) # Cholesky factorisation on sparse symmetric matrix

w = zeros(J+1)

for n ∈ 1:N
    w[2:end-1] = u[2:end-1] - .5σ * (u[3:end] - u[1:end-2])
    w[1] = u0(x[1] - V * n * Δt)
    w[end] = w[1]
    u = F\w
end

uexact = u0.(x .- V * tf)
err = u - uexact
u_plot = plot(x, u, xlabel="x", ylabel="u", color=:red, label="Finite differences", lw=6)
plot!(u_plot, x, uexact, xlabel="x", ylabel="u", color=:black, label="Solution", lw=2)
err_plot = plot(x, err, xlabel="x", ylabel="Error", legend=false)
display(plot(u_plot, err_plot, layout=(2, 1), size=(700, 700)))
println("Δx: ", Δx, "\t Δt:", Δt, "\t max error: ", maximum(abs.(err)))
```