clear all;
close all;
function U = RoomTemperature(ot,dt,ht,n)
% ROOMTEMPERATURE calcule la température ambiante dans une chambre
% U = RoomTemperature(ot, dt, ht, n); prend la température extérieure ot,
% la température de la porte dt et la température du chauffage ht et la
% nombre de points de grille n, puis calcule la température ambiante de la
% chambre
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = (((X>-0.5) & (X<0.5) & (Y>-0.5) & (Y<0.5)) |...
    ((X>-1) & (X<-0.5) & (Y>-0.2) & (Y<0.2))|...
    ((X>0.5) & (X<1) & (Y>-0.2) & (Y<0.2)));


H = ((X>-0.2) & (X<0.2) & (Y>-0.2) & (Y<0.2));      % chauffage au centre

k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);% Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              % CL Neumann sur le mur gauche
        if (Y(i,j)>-0.2 && Y(i,j)<0.2) 
          door = [door ind];
        else
        A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i,j+1) == 0 && j<n-1         % CL Neumann sur le mur droit
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0      % On garde la CL Dirichlet pour la fenetre
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
      end
      if G(i-1,j) == 0
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
      end 
    end
  end
end
h = 2/(n-1);
A = -A/h^2;                         % division par h^2
b = zeros(length(k),1);
window = G(G(:,end-1)>0,end-1);     % indices correspondant a la fenetre


heat = G(H);                        % indices correspondant au radiateur faire attention qu'il se trouve bien dans la pièce

b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(heat) = -ht;                      % et du chauffage
b(door) = -1/h^2*dt;
I = eye(size(A));

alpha = 0.5; nu = 1;
deltat = alpha/(2*nu)*(h^2);


u = ot*ones(length(k),1); %initialisation de la température de la chambre
for l = 1:10000
  uold = u;
  %u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); %implicite
  u = (eye(size(A))+nu*deltat*A)*uold - nu*deltat*b; %explicite
  
  U = G;
  U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
  if mod(l,500) == 0   %on affiche toutes les 500 itérations
    mesh(X,Y,U);                        % correspondant a la grille 
    axis('ij');                         % on dessine la solution sur la grille
    pause(0.1);
  end
end                    % on dessine la solution sur la grille
disp(mean(mean(U((((X>-0.5) & (X<0.5) & (Y>-0.5) & (Y<0.5)) |...    % Pour avoir la temperature Moyenne de la chambre
    ((X>-1) & (X<-0.5) & (Y>-0.2) & (Y<0.2))|...                      
    ((X>0.5) & (X<1) & (Y>-0.2) & (Y<0.2)))))));
end
U=RoomTemperature(-10,8,170,40);