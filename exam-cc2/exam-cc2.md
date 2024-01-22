![PNS](http://caillau.perso.math.cnrs.fr/logo-pns.png)
## MAM4
# Équations aux dérivées partielles 1
# 2023-24

# Exam CC no. 2

**Durée 2H. Documents non autorisés. Tous les exercices sont indépendants. Le barème indicatif est précisé pour chaque exercice.**

## Exercice 1 (5 points)

### 1.1
Soit $w : \mathbf{R}^n \to \mathbf{R}^n$ une fonction suffisamment régulière, et $\Omega$ un ouvert de $\mathbf{R}^n$ à bord $\Gamma := \partial\Omega$ également régulier. Démontrer la *formule de la divergence* à l'aide de la formule de Stockes :

$$ \int_\Omega \nabla \cdot w\ \mathrm{d}x = \int_\Gamma w \cdot n\ \mathrm{d}\sigma $$

où $\nabla = (\partial/\partial x_1,\dots,\partial/\partial x_n)$ et où $n$ désigne la normale extérieure au bord $\Gamma$.

**Réponse.**

$$ \sum_{i=1}^n \int_\Omega \frac{\partial w_i}{\partial x_i}\ \mathrm{d}x
 = \sum_{i=1}^n \int_\Gamma w_i n_i\ \mathrm{d}\sigma $$

### 1.2
On considère l'équation des ondes en dimension un d'espace avec conditions aux limites périodiques :

$$ \frac{\partial^2 u}{\partial t^2}(x,t) - \frac{\partial^2 u}{\partial x^2}(x,t) = 0,\quad x \in \mathbf{R},\quad t > 0, $$

$$ u(x+1,t) = u(x,t), \quad x \in \mathbf{R},\quad t > 0, $$

$$ u(x,0) = u_0(x),\quad \frac{\partial u}{\partial t}(x,0) = u_1(x). $$

Montrer que, pour $x$ dans $\mathbf{R}$ et $t > 0$, on a

$$ \frac{\partial u}{\partial x}(x+1,t) = \frac{\partial u}{\partial x}(x,t), $$

et

$$ \frac{\partial u}{\partial t}(x+1,t) = \frac{\partial u}{\partial t}(x,t). $$

En déduire que, pour toute solution suffisamment régulière, *l'énergie* ci-dessous est conservée au cours du temps :

$$ E(t) := \frac{1}{2} \int_0^1 \left( \frac{\partial u}{\partial x} \right)^2 + \left( \frac{\partial u}{\partial t} \right)^2 \ \mathrm{d}x. $$

**Réponse.** On obtient les relations voulues en dérivant soit par rapport à $x$, soit par rapport à $t$ l'égalité $u(x+1,t) = u(x,t)$. Pour une solution $u$ suffisamment régulière, on peut ensuite dériver par rapport au temps sous le signe somme

$$ E'(t) = \int_0^1 \left\( \frac{\partial u}{\partial x} \frac{\partial^2 u}{\partial t\partial x} + \frac{\partial u}{\partial t} \frac{\partial^2 u}{\partial t^2}\right\)\ \mathrm{d}x, $$

puis intégrer par parties selon

$$ \int_0^1 \frac{\partial u}{\partial x} \frac{\partial^2 u}{\partial t\partial x}\ \mathrm{d}x = \left[ \frac{\partial u}{\partial x} \frac{\partial u}{\partial t} \right]_0^1 - \int_0^1 \frac{\partial^2 u}{\partial x^2} \frac{\partial u}{\partial t}\ \mathrm{d}x. $$

En utilisant les relations précédentes, on voit que le terme intégré est nul, et les termes restant s'annulent parce que $u$ est solution.

## Exercice 2 (8 points)
On considère le problème aux limites à condition de Neumann hétérogène sur $\Omega$, ouvert borné connexe à bord $\Gamma := \partial\Omega$ régulier de $\mathbf{R}^n$ : trouver $u$ dans $\mathscr{C}^2(\overline{\Omega})$ telle que (*solution forte*)

$$ -\Delta u(x) + \gamma u(x) = f(x),\quad x \in \Omega, $$

$$ \frac{\partial u}{\partial n}(\sigma) = g(\sigma),\quad \sigma \in \Gamma, $$

où $f$ appartient à $\mathscr{C}^0(\overline{\Omega})$ et $g$ à  $\mathscr{C}^0(\Gamma)$, et où $\gamma$ est un réel strictement positif.

### 2.1
On pose $H := H^1(\Omega)$. En la testant contre une fonction arbitraire $v$ de $H$, montrer que toute solution forte est nécessairement solution d'un nouveau problème (*solution faible*) de la forme

$$ a(u,v) = \varphi \cdot v,\quad v \in H, $$

où $a$ et $\varphi$ sont des formes bilinéaire et linéaire, respectivement, que l'on précisera.

**Réponse.** On utilise Stockes (licite car $u$ est de classe $\mathscr{C}^2$ et $v$ est dans $H^1$) pour écrire 

$$ -\int_\Omega \Delta u v\ \mathrm{d}x = -\int_\Gamma \frac{\partial u}{\partial n} v\ \mathrm{d}\sigma + \int_\Omega \nabla u \cdot \nabla v\ \mathrm{d}x. $$

Comme $\partial u/\partial n = g$ sur $\Gamma$, on obtient la forme voulue avec

$$ a(u,v) := \int_\Omega \nabla u \cdot \nabla v\ \mathrm{d}x + \gamma \int_\Omega u v\ \mathrm{d}x \text{ et } \varphi(v) := \int_\Omega f v\ \mathrm{d}x + \int_\Gamma g v\ \mathrm{d}\sigma. $$

### 2.2
Montrer que $a$ est continue et coercive.

**Réponse.** Pour $u$ et $v$ dans $H$ on a

$$ |a(u,v)| \leq (1+\gamma) ||u||_{H^1} \|v\| {H^1} $$

d'où la continuité, et

$$ a(u,u) \geq \min(1,\gamma) \|u\|^2_{H^1} $$

d'où la coercivité.

### 2.3
Montrer que $\varphi$ est continue.

**Nota bene**. On rappelle que l'opérateur de trace $v \mapsto v_{|\Gamma}$ de $\mathscr{C}^1(\overline{\Omega})$ dans $L^2(\Gamma)$, qui associe à une fonction sa restriction sur le bord, se prolonge continûment sur tout $H^1(\Omega)$. 

**Réponse.** Par Cauchy-Schwarz

$$ |\varphi(v)| \leq \|f\|_{L^2} \|v\|_{L^2} + \|g\|_{L^2} \|v_{|\Gamma}\|_{L^2}, $$

et, par continuité de la trace, il existe une constante $C$ telle que $\|v_{|\Gamma}\|_{L^2} \leq C \|v\|_{H^1}$, d'où la continuité de $\varphi$.

### 2.4
En déduire l'existence et l'unicité de solution faible, et montrer que cette solution faible est également solution d'un problème variationnel que l'on précisera.

**Réponse.** L'existence et l'unicité de solution faible découlent du théorème de Lax-Milgram. Comme $a$ est symétrique, le même résultat assure en outre que la solution faible est également solution du problème variationnel

$$ \frac{1}{2}a(u,u) - \varphi \cdot u \to \min,\quad u \in H. $$

### 2.5
On suppose que la solution faible est de classe $\mathscr{C}^2$ sur $\overline{\Omega}$. Montrer qu'elle est aussi solution forte.

**Réponse.** Pour $v$ arbitraire dans $H$, sous l'hypothèse que $u$ est de classe $\mathscr{C}^2$ on peut utiliser la formule de Stockes sur la formulation faible pour écrire

$$ \int_\Omega \nabla u \cdot \nabla v\ \mathrm{d}x = \int_\Gamma \frac{\partial u}{\partial n} v\ \mathrm{d}\sigma - \int_\Omega \Delta u v\ \mathrm{d}x, $$

et montrer que

$$ \int_\Omega (-\Delta u + \gamma u - f) v\ \mathrm{d}x + \int_\Gamma (\frac{\partial u}{\partial n} - g)v\ \mathrm{d}\sigma = 0. $$

En prenant $v$ nulle au bord et en utilisant le fait que $H^1_0$ est dense dans $L^2$, on obtient $-\Delta u + \gamma u = f$ dans $\Omega$ (partout, $u$ étant $\mathscr{C}^2$). Par conséquent, pour tout $v$ dans $H$,

$$ \int_\Gamma (\frac{\partial u}{\partial n} - g)v\ \mathrm{d}\sigma = 0, $$

d'où l'on déduit que $\partial u/\partial n = g$ sur $\Gamma$ (densité de l'image de l'application trace dans $L^2(\Gamma)$).

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

**Réponse.** On fait un DL en $(x_j,t_n)$ d'une solution $u$ suffisamment régulière, à l'ordre deux en temps et trois en espace pour obtenir que le schéma est consistant d'ordre un en temps et deux en espace :

$$ E_j^n = \frac{1}{2}\frac{\partial^2 u}{\partial t^2}(x_j,t_n) \Delta t + O(\Delta t^2) + \frac{V}{6}\frac{\partial^3 u}{\partial x^3}(x_j,t_n) \Delta x^2 + O(\Delta x^3). $$

### 3.2
Montrer que ce schéma est inconditionnellement instable au sens $L^2$.

**Réponse.** On injecte le mode $u = A^n(k) e^{2i\pi kj\Delta x}$ dans le schéma pour obtenir, après simplification ($\sigma := V\Delta t/\Delta x$)),

$$ |A(k)|^2 = 1 + \sigma^2\sin^2(2\pi k\Delta x) \geq 1. $$

### 3.3
Quel schéma obtient-on si on utilise l'équation équivalente pour améliorer le schéma précédent ?

**Réponse.** En retrenchant le terme principal (d'ordre un) de l'erreur,

$$ frac{1}{2}\frac{\partial^2 u}{\partial t^2}(x_j,t_n) \Delta t = \frac{V^2}{2}\frac{\partial^2 u}{\partial x^2}(x_j,t_n) \Delta t \simeq \frac{V^2 \Delta t}{2}\frac{u(x_{j+1},t_n)-2u(x_j,t_n)+u(x_{j-1},t_n}{\Delta x^2}, $$

on obtient le schéma de Lax-Wendroff explicite

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + \frac{V}{2} \frac{u_{j+1}^n-u_{j-1}^n}{\Delta x} - \frac{V^2 \Delta t}{2}\frac{u_{j+1}^n-2u_j^n+u_{j-1}^n}{\Delta x^2} = 0. $$

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

**Réponse.** Le schéma de Lax-Wendroff implicite (voir TP Advection),

$$ \frac{u_j^{n+1}-u_j^n}{\Delta t} + \frac{V}{2} \frac{u_{j+1}^{n+1}-u_{j-1}^{n+1}}{\Delta x}
 - \frac{V^2 \Delta t}{2}\frac{u_{j+1}^{n+1}-2u_j^{n+1}+u_{j-1}^{n+1}}{\Delta x^2} = 0 $$

où $\sigma := V\Delta t/\Delta x$.

### 3.5
Toujours dans ce même code, expliquer le rôle de la ligne `u = F\w`. En particulier, que désigne `F` ? Comment peut-on améliorer la performance de  la boucle `for n...`  ?

**Réponse.** On stocke dans $F$ la factorisation de Cholesky de la matrice $A$. La matrice $A$ étant indépendante de $n$, on doit extraire cette factorisation de la boucle pour ne la faire qu'une seule fois.