% Simulation du probleme de Poisson
% Probleme de Dirichlet non-homogene
% Avec un terme source localise (de type chauffage)
% Prise en compte des CL Dirichlet (ajout d'une fenetre + porte)
ht = 100; % temperature du chauffage
dt = 15;  % temperature de la porte d'entree
wt = -10; % temperature de la fenetre
n = 30;   % nombre de points de discretisation

x = linspace(-1,1,n);                      
[X,Y] = meshgrid(x,x);                  % grille bi-dimensionnelle 
G = ((X>-1) & (X<1) & (Y>-1) & (Y<1));  % domaine: on prend le carre en entier
H = ((X>-0.6) & (X<0.5) & (Y>-0.1) & (Y<0.1)); % position du chauffage

k = find(G);                          % les indices des G non nuls
G = zeros(size(G));                   % conversion logique - reel
G(k) = (1:length(k))';                % suite d'indices de la matrice A
spy(G)                                % on dessine la grille
heat = G(H);                          % indices correspondant a H

A = delsq(G);                         % matrice du Laplacien a l'interieur
door = [];
window = [];
for i = 2:n-1                        
  for j = 2:n-1                      
    if G(i,j) ~= 0 
      if G(i+1,j) == 0                % la porte se situe sur le mur du haut
        if (X(i,j)> -0.9 && X(i,j)< -0.5) 
          door = [door G(i,j)];       % indices correspondant a la porte 
        end
      end
      if G(i,j+1) == 0                % la fenetre se situe sur le mur a droite
        if (Y(i,j)>0.5 && Y(i,j)<0.9) 
          window = [window G(i,j)];   % indices correspodant a la fenetre 
        end
      end
    end
  end
end

h = 2/(n-1);                          % pas du maillage
A = -A/h^2;                           % division par h^2
b = zeros(length(k),1);               % initialisation du second membre
b(heat) = -ht;                        % prise en compte du chauffage
b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*wt;                % CL Dirichet: fenetre
u = A\b;                              % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice 
                                      % correspondant a la grille
mesh(X,Y,U);                          % on dessine la solution sur la grille
axis('ij');