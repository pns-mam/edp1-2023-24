clear all, close all, clc
ot=-10;
dt=15;
ht=20;                                        %situation en hiver
n=40;
%ot=30;                                       %situation en ete
%dt=20;
%ht=20;                                              
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-0.8) & (X<0.6) & (Y>-0.4) & (Y<0.6)) | ...  % geometrie de la chambre
   ((X>-0.2) & (X<0.6) & (Y>0.6) & (Y<0.9)) | ... 
   ((X>-0.4) & (X<0.6) & (Y>-0.7) & (Y<-0.4));
%H = ((X>-0.2) & (X<0.2) & (Y>-0.1) & (Y<0.1));       % position du chauffage au centre
H = ((X>0.4) & (X<0.6) & (Y>0.2) & (Y<0.6));      % chauffage pres le mur
%H = ((X>0) & (X<0.4) & (Y>-0.7) & (Y<-0.5));      % chauffage pres le fenetre
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              % CL Neumann sur le mur gauche
          if (Y(i,j)>0.1 && X(i,j)<0.5) 
           door = [door ind];        % On garde la CL Dirichlet pour la porte 
          else
           A(ind,ind) = A(ind,ind)-1;
          end
      end
      if G(i,j+1) == 0              % CL Neumann sur le mur droit
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0 && i<n-1     % On garde la CL Dirichlet pour la fenetre
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
window = G(end-1,G(end-1,:)>0);     % indices correspondant a la fenetre
heat = G(H);                        % indices correspondant au radiateur
b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(heat) = -ht;                      % et du chauffage
b(door) = -1/h^2*dt;
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');                         % on dessine la solution sur la grille