function [U,G1] = RoomTemperature2(ot,dt,ht,n,titre)
% ROOMTEMPERATURE calcule la température ambiante dans une chambre
% U = RoomTemperature(ot, dt, ht, n); prend la température extérieure ot,
% la température de la porte dt et la température du chauffage ht et la
% nombre de points de grille n, puis calcule la température ambiante de la
% chambre

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G1 = ((X>-0.7) & (X<0.7) & (Y>0) & (Y<0.9)) | ((X>-0.5) & (X<0.5) & (Y>-0.8) & (Y<0)) ;

H = ((X>-0.2) & (X<0.2) & (Y>-0.5) & (Y<-0.2));
k = find(G1);                        % les indices des G non nuls
G = zeros(size(G1));                 % conversion logique - reel 
G(k) = (1:length(k))';              % suite d'indices de la matrice A  

A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
window = [];
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i+1,j) == 0              % mur du haut isolant donc Neumann
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i-1,j) == 0              % mur du bas
        if (X(i,j)>-0.25 && X(i,j) < 0.25)
          window = [window ind]; %fenetre
        else
          A(ind,ind) = A(ind,ind)-1; %sinon isolant
        end
        
      end
      if G(i,j-1) == 0  %mur a gauche
        if (Y(i,j) >0.3 && Y(i,j) < 0.8)
          door = [door ind];
        else
          A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
        end
        
      end
      if G(i,j+1) == 0 %mur a droite
          A(ind,ind) = A(ind,ind)-1;
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
U(G>0) = full(u(G(G>0))); 
figure          % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');  
title(titre)
  