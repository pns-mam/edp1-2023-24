% Simulation du probleme de Poisson dans le cas du probleme de Dirichlet homogene
n = 30;  % nombre de points de discretisation

x = linspace(-1,1,n);                      
[X,Y] = meshgrid(x,x);                  % grille bi-dimensionnelle 
G = ((X>-1) & (X<1) & (Y>-1) & (Y<1));  % on prend le carre en entier      
k = find(G);                          % On trouve les indices des G non nuls
G = zeros(size(G));                   % Conversion logique - reel
G(k) = (1:length(k))';                % Suite d'indices de la matrice A
spy(G)
pause

A = delsq(G);                       % Matrice du -Laplacien sans le 1/h^2
h = 2/(n-1);                          % pas du maillage
A = -A/h^2;                           % Division par h^2
b = -ones(length(k),1);               % initialisation du second membre

u = A\b;                              % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice 
                                      % correspondant a la grille
mesh(X,Y,U);                          % on dessine la solution sur la grille
axis('ij');