% Simulation du probleme de Poisson
% Probleme de Dirichlet non-homogene
% Avec un terme source localise (de type chauffage)
% Prise en compte des CL Dirichlet (ajout d'une fenetre + porte)

ht = 900; % temperature du chauffage
dt = 15;  % temperature de la porte d'entree
ot = -10; % temperature de la fenetre
n = 30;   % nombre de points de discretisation

x = linspace(-3,3,n);
[X,Y] = meshgrid(x,x);                % grille bi-dimensionnelle
G = (((X>-3) & (X<0) & (Y>-1) & (Y<2)) | ((X>0) & (X<2) & (Y>-1) & (Y<1))); % Dimensions de la Chambre
H = (((X>-3) & (X<-2.5) & (Y>-1) & (Y<1)) | ((X>1.5) & (X<2) & (Y>0) & (Y<1))); % position des chauffages
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
window = [];
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              % CL Neumann sur le mur gauche
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i,j+1) == 0              % CL Neumann sur le mur droit
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0 && !((X(i,j)>-3 && X(i,j)<-1)|(X(i,j)>0.5 && X(i,j)<1.5))
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas
      end
      if G(i+1,j) == 0 &&((X(i,j)>-3 && X(i,j)<-1)|(X(i,j)>0.5 && X(i,j)<1.5)) % On garde la CL Dirichlet pour les fenetres (situées sur le mur du bas)
        window = [window ind];
      end
      if G(i-1,j) == 0
        if (X(i,j)>-2.5 && X(i,j)<-1.7)
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        else
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
        end
      end 
    end
  end
end
h = 2/(n-1);
A = -A/h^2;                         % division par h^2
b = zeros(length(k),1);
heat = G(H);                        % indices correspondant au radiateur
b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(heat) = -ht;                      % et du chauffage
b(door) = -1/h^2*dt;
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');                         % on dessine la solution sur la grille

min(min(U))
max(max(U))
m=0;
for i=1:length(k)
  m+=u(i);
end
m/length(k)
nb=0;
for i=1:length(k)
  if 18<= u(i) && u(i)<=23
    nb+=1;
  end
end
nb/length(k)

