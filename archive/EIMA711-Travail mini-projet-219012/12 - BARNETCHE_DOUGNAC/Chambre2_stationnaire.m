%{
MINI PROJET - EDP1 MAM4
Modélisation de la température dans une pièce
DOUGNAC Jade / BARNETCHE Carine
%}

% PARTIE 1 : CHAMBRE 2 - MODELE STATIONNAIRE


clear all, close all, clc
function U = RoomTemperature(ot,dt,ht,n)
% ROOMTEMPERATURE calcule la température ambiante dans une chambre
% U = RoomTemperature(ot, dt, ht, n); prend la température extérieure ot,
% la température de la porte dt et la température du chauffage ht et la
% nombre de points de grille n, puis calcule la température ambiante de la
% chambre
  x = linspace(-1,1,n);                            
  [X,Y] = meshgrid(x,x);                           % grille bi-dimensionnelle 
  G = ((X>-1) & (X<0.3) & (Y>-1) & (Y<0.8)) | ...  % geometrie de la chambre
      ((X>0.3) & (X<0.9) & (Y>0) & (Y<0.8));

  H = ((X>-0.4) & (X<0) & (Y>-0.6) & (Y<-0.3));     % chauffage

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
        if G(i,j-1) == 0            
          if (Y(i,j+1)>-0.5 && Y(i,j+1)<0.1) 
            window = [window ind];     % On garde la CL Dirichlet pour la fenêtre
          else
             A(ind,ind) = A(ind,ind)-1;% CL Neumann sur le mur gauche
          end
        end
        if G(i,j+1) == 0       
          A(ind,ind) = A(ind,ind)-1;   % CL Neumann sur le mur droit
        end
        if G(i+1,j) == 0      
          A(ind,ind) = A(ind,ind)-1;   % CL Neumann sur le mur du bas 
        end
        if G(i-1,j) == 0 
          if (X(i,j)>0.4 && X(i,j)<0.7) 
            door = [door ind];         % On garde la CL Dirichlet pour la porte 
         else
            A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
          end
        end 
      end
    end
  end

  h = 2/(n-1);                        % pas d'espace
  A = -A/h^2;                         % division par h^2
  b = zeros(length(k),1);

  heat = G(H);                        % indices correspondant au radiateur

  b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
  b(door) = -1/h^2*dt;
  b(heat) = -ht;                      % et du chauffage

  u = A\b;                            % solution du systeme sous forme de vecteur
  U = G;
  U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
  mesh(X,Y,U);                        % correspondant a la grille 
  axis('ij');                         % on dessine la solution sur la grille
end

U = RoomTemperature(-10,15,300,40); 

% On affiche la température moyenne dans la pièce
% on fait bien attention à calculer les valeurs à l'intérieur de la chambre
x = linspace(-1,1,40);                            
[X,Y] = meshgrid(x,x);
temperature_moyenne = mean(mean(U((((X>-1) & (X<0.3) & (Y>-1) & (Y<0.8)) | ...  % geometrie de la chambre
    ((X>0.3) & (X<0.9) & (Y>0) & (Y<0.8))))));
disp(temperature_moyenne);