
tic
ot=-10;
dt=15;
ht=500;
n=40;
% Nous allons calculer la température ambiante dans la chambre 2
% a l'aide de l'equation de la chaleur en mode stationnaire.
% ot est la température extérieure et donc celle de la fenêtre,
% dt est la température de la porte et ht est la température du chauffage et
% n est le nombre de points de grille.

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
   
G = ((X>-0.3) & (X<0.8) & (Y>-0.9) & (Y<0.8)) | ...  % geometrie de la chambre
   ((X>-0.8) & (X<-0.3) & (Y>-0.9) & (Y<-0.5)) | ...
   ((X>-0.8) & (X<0.2) & (Y>-0.9) & (Y<-0.5)) | ...
   ((X<0.95) & (X>0.8) & (Y>-0.3) & (Y<0.3));

H = ((X<0.75) & (X>0.5) & (Y>-0.3) & (Y<0.3));       % position du chauffage
                                                    % en dessous de la fenetre
%H = ((X>-0.1) & (X<0.5) & (Y>0.5) & (Y<0.75));      % position du chauffage
                                                     % sur le mur du fond                                                    
%H = ((X>-0.25) & (X<0) & (Y>-0.3) & (Y<0.3));      % position du chauffage
                                                     % sur le mur de gauche                                                     
                                                     
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
      if G(i+1,j) == 0      
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
      end
      if G(i,j+1) == 0                     % la fenetre se situe sur le mur a droite
          if (Y(i,j)>-0.3 && Y(i,j)<0.3) 
            window = [window G(i,j)];   % indices correspodant a la fenetre 
          else
            A(ind,ind) = A(ind,ind)-1;
          end
      end
      if G(i-1,j) == 0
        if (Y(i,j)>-0.9 && Y(i,j)<-0.5) 
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
title('Chambre 2 avec chauffage au niveau de la fenetre');
toc