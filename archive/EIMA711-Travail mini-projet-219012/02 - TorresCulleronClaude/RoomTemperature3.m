function [U,G1] = RoomTemperature3(ot,dt,ht,n,titre)
% ROOMTEMPERATURE calcule la température ambiante dans une chambre
% U = RoomTemperature(ot, dt, ht, n); prend la température extérieure ot,
% la température de la porte dt et la température du chauffage ht et la
% nombre de points de grille n, puis calcule la température ambiante de la
% chambre

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G1 = ((X>-1) & (X<1) & (Y>-0.33) & (Y<0.8)) | ...  % geometrie de la chambre
   ((X>-1) & (X<-0.5) & (Y>-1) & (Y<-0.33)) | ... 
   ((X>-0.6) & (X<0) & (Y>0.8) & (Y<1));

H = ((X>-0.6) & (X<-0.2) & (Y>0.5) & (Y<0.8)); %position améliorée du chauffage
k = find(G1);                        % les indices des G non nuls
G = zeros(size(G1));                 % conversion logique - reel
disp('Val non nulle');
disp(length(k)); 
 
G(k) = (1:length(k))';              % suite d'indices de la matrice A 
 
A = delsq(G);                       % Laplacien a l'interieur du domaine G
disp(size(A))
door = [];                          % indices correspondant a la porte
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
      if G(i+1,j) == 0 && i<n-1     % On garde la CL Dirichlet pour la fenetre
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
      end
      if G(i-1,j) == 0
        if (X(i,j)>-0.9 && X(i,j)<-0.7) 
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

window = G(end-1,G(end-1,:)>0);     % indices correspondant a la fenetre

heat = G(H);                        % indices correspondant au radiateur
b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(heat) = -ht;                      % et du chauffage
b(door) = -1/h^2*dt;
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;

U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
figure
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');                         % on dessine la solution sur la grille
title(titre)