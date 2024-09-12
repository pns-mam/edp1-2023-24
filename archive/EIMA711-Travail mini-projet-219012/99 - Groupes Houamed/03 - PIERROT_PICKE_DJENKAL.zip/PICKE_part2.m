% Simulation du probleme de Poisson
% Probleme de Dirichlet non-homogene
% Avec un terme source localise (de type chauffage)
% Prise en compte des CL Dirichlet (ajout d'une fenetre + porte)
ht = -1800; % temperature du chauffage
dt = 30;  % temperature de la porte d'entree
wt = 30; % temperature de la fenetre
n = 30;   % nombre de points de discretisation

tic
x = linspace(-3,3,n);
[X,Y] = meshgrid(x,x);                % grille bi-dimensionnelle
G = (((X>-3) & (X<0) & (Y>-1) & (Y<2)) | ((X>0.1) & (X<2) & (Y>-1) & (Y<1)));
H = (((X>-3) & (X<-2.5) & (Y>-1) & (Y<1)) | ((X>1.5) & (X<2) & (Y>0) & (Y<1))); % position du chauffage


k = find(G);                          % les indices des G non nuls
G = zeros(size(G));                   % conversion logique - reel
G(k) = (1:length(k))';                % suite d'indices de la matrice A
spy(G)                                % on dessine la grille
heat = G(H);                         % indices correspondant a H

A = delsq(G);                       % matrice du Laplacien a l'interieur
door = [];
window = [];
for i = 2:n-1
    for j = 2:n-1
        if G(i,j) ~= 0
            if G(i-1,j) == 0                % la porte se situe sur le mur du haut
                if (X(i,j)> -2.5 && X(i,j)< -1.7)
                    door = [door G(i,j)];       % indices correspondant a la porte
                end
            end
            if G(i+1,j) == 0                % la fenetre se situe sur le mur a droite
                if ((X(i,j)>-3 && X(i,j)<-1)||(X(i,j)>0.5 && X(i,j)<1.5))
                    window = [window G(i,j)];   % indices correspodant a la fenetre
                end
            end
        end
    end
end

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
      if G(i+1,j) == 0 && !((X(i,j)>-3 && X(i,j)<-1)|(X(i,j)>0.5 && X(i,j)<1.5))    % On garde la CL Dirichlet pour la fenetre
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas
      end
      if G(i-1,j) == 0
        if !(X(i,j)>-2.5 && X(i,j)<-1.7) 
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
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

alpha = 0.25; nu =1;
deltat = alpha/(nu)*h^2;
u = 30*ones(length(k),1);     % solution initiale

for l = 1:400           % boucle en temps 
    uold = u;
    %u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
    %u = (eye(size(A))+nu*deltat*A)*uold - nu*deltat*b; % explicite
    u=(eye(size(A))-nu*deltat*A)\((eye(size(A))+nu*deltat*A)*uold-nu*deltat*b);% Crank Nicolson
    u(door) = dt;                  % ET ICI
    u(window) = wt;                % ET ICI
    U = G;
    U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice correspondant à la grille
    %mesh(X,Y,U);                          % on dessine la solution sur la grille
    %axis('ij');
    %pause(0.005);
end

mesh(X,Y,U);                          % on dessine la solution sur la grille
axis('ij');
toc

min=u(1);
max=u(1);
m=0;
for i=1:length(k)
  m+=u(i);
  if (u(i)<min)
    min=u(i);
   end
   if (u(i)>max)
     max=u(i);
    end
end
min
max

moyenne=m/length(k)
nb=0;
for i=1:length(k)
  if 18<= u(i) && u(i)<=23
    nb+=1;
  end
end
p=nb/length(k)

"Fini"