![PNS](http://caillau.perso.math.cnrs.fr/logo-pns.png)
## MAM4
# Équations aux dérivées partielles 1
# 2023-24

# Exam CC no. 2

**Durée 2H. Documents non autorisés. Tous les exercices sont indépendants. Le barème indicatif est précisé pour chaque exercice.**

## Exercice 1 (5 points)

### 1.1
Soit $w : \mathbf{R}^n \to \mathbf{R}^n$ une fonction suffisamment régulière, et $\Omega$ un ouvert de $\mathbf{R}^n$ à bord $\Gamma := \partial\Omega$ également régulier. Démontrer la *formule de la divergence* à l'aide de la formule de Green :

$$ \int_\Omega \nabla \cdot w\ \mathrm{d}x = \int_\Gamma w \cdot n\ \mathrm{d}\sigma $$

où $\nabla = (\partial/\partial x_1,\dots,\partial/\partial x_n)$ et où $n$ désigne la normale extérieure au bord $\Gamma$.


### 1.2
On considère l'équation des ondes en dimension un d'espace avec conditions aux limites périodiques :

$$ \frac{\partial^2 u}{\partial t^2}(x,t) - \frac{\partial^2 u}{\partial x^2}(x,t) = 0,\quad x \in \mathbf{R},\quad t > 0, $$

$$ u(x+1,t) = u(x,t), \quad x \in \mathbf{R},\quad t > 0, $$

$$ u(x,0) = u_0(x),\quad \frac{\partial u}{\partial t}(x,0) = u_1(x). $$

Montrer que, pour $x$ dans $\mathbf{R}$ et $t > 0$, on a

$$ \frac{\partial u}{\partial x}(x+1,t) = \frac{\partial u}{\partial x}(x,t), $$

et

$$ \frac{\partial u}{\partial t}(x+1,t) = \frac{\partial u}{\partial t}(x,t). $$

En déduire que *l'énergie* ci-dessous est conservée au cours du temps :

$$ E(t) := \frac{1}{2} \int_0^1 \left( \frac{\partial u}{\partial x} \right)^2 + \left( \frac{\partial u}{\partial t} \right)^2 \ \mathrm{d}x. $$

## Exercice 2 (8 points)
On considère le problème aux limites à condition de Neumann hétérogène sur $\Omega$, ouvert borné connexe à bord $\Gamma := \partial\Omega$ régulier de $\mathbf{R}^n$ : trouver $u$ dans $\mathscr{C}^2(\overline{\Omega})$ telle que (*solution forte*)

$$ -\Delta u(x) + \gamma u(x) = f(x),\quad x \in \Omega, $$

$$ \frac{\partial u}{\partial n}(\sigma) = g(\sigma),\quad \sigma \in \Gamma, $$

où $f$ appartient à $\mathscr{C}^0(\overline{\Omega})$ et $g$ à  $\mathscr{C}^0(\Gamma)$, et où $\gamma$ est un réel strictement positif.

### 2.1
On pose $H := H^1(\Omega)$. En la testant contre une fonction arbitraire $v$ de $H$, montrer que toute solution forte est nécessairement solution d'un nouveau problème (*solution faible*) de la forme

$$ a(u,v) = \varphi \cdot v,\quad v \in H, $$

où $a$ et $\varphi$ sont des formes bilinéaire et linéaire, respectivement, que l'on précisera.

### 2.2
Montrer que $a$ est continue et coercive.

### 2.3
Montrer que $\varphi$ est continue.

**Nota bene**. On rappelle que l'opérateur de trace $v \mapsto v_{|\Gamma}$ de $\mathscr{C}^1(\overline{\Omega})$ dans $L^2(\Gamma)$, qui associe à une fonction sa restriction sur le bord, se prolonge continûment sur tout $H^1(\Omega)$. 

### 2.4
En déduire l'existence et l'unicité de solution faible, et montrer que cette solution faible est également solution d'un problème variationnel que l'on précisera.

### 2.5
On suppose que la solution faible est de classe $\mathscr{C}^2$ sur $\overline{\Omega}$. Montrer qu'elle est aussi solution forte.

## Exercice 3 (7 points)

On considère l'équation d'advection en dimension un d'espace avec conditions aux limites périodiques :

$$ \frac{\partial u}{\partial t}(x,t) + V \frac{\partial u}{\partial x}(x,t) = 0,\quad x \in ]0,1[,\quad t > 0, $$

$$ u(0,t) = u(1,t),\quad t > 0, $$

$$ u(x, 0) = u_0(x),\quad x \in ]0,1[. $$

On suppose $V > 0$.

### 3.1
On discrétise le problème à l'aide du schéma explicite centré ci-dessous :

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + \frac{V}{2} \frac{u_{j+1}^n-u_{j-1}^n}{\Delta x} = 0. $$

Montrer que ce schéma est consistant et préciser les ordres d'approximation en espace et en temps.

### 3.2
Montrer que ce schéma est inconditionnellement instable au sens $L^2$.

### 3.3
Quel schéma obtient-on si on utilise l'équation équivalente pour améliorer le schéma précédent ?

### 3.4
Dans le code ci-dessous, quel schéma numérique a-t-on implémenté ?

```julia
u = u0(x)
A = spdiagm(-1 => -σ^2/2*ones(Nx), 0 => (1+σ^2)*ones(Nx+1), 1 => -σ^2/2*ones(Nx))
A[1, end] = -σ^2/2
A[end, 1] = -σ^2/2 
w = zeros(Nx+1)

for n ∈ 1:Nt
    w[2:end-1] = u[2:end-1] - σ/2 * (u[3:end]-u[1:end-2])
    w[1] = u0(x[1] - V * n * Δt)
    w[end] = w[1]
    F = cholesky(A)
    u = F\w
end
```

### 3.5
Toujours dans ce même code, expliquer le rôle de la ligne `u = F\w`. En particulier, que désigne `F` ? Comment peut-on améliorer la performance de  la boucle `for n...`  ?