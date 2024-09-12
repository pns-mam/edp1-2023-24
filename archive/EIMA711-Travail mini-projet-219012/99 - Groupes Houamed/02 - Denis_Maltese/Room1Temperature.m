function U = Room1Temperature(ot,dt,ht,n)
% ROOMTEMPERATURE calcule la température ambiante dans la chambre 1
% ot : température extérieure ot
% dt : température de la porte  
% ht : température du chauffage ht 
% n : nombre de points de grille n

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-1) & (X<1) & (Y>0) & (Y<1)) | ((X>-1) & (X<0.2) & (Y>-1) & (Y<1));
%géométrie de la chambre
H = ((X>-0.7) & (X<-0.1) & (Y>-0.8) & (Y<-0.5));      % position du chauffage
%H=((X>-0.7) & (X<0.1) & (Y>0.2) & (Y<0.5));          %chauffage au centre
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
window2 = [];                       % indices correspondant a la fenetre du bas
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
      if G(i+1,j) == 0 
        if (X(i,j)>-1 && X(i,j)<-0.2)   % On garde la CL Dirichlet pour la fenetre 2
          window2 = [window2 G(i,j)];
        else   
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
        end
      end
      if G(i-1,j) == 0 && i>2       % On garde la CL Dirichlet pour la fenetre 1
        if (X(i,j)>0.4 && X(i,j)<0.8) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
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
window1 = G(2,G(2,:)>0);     % indices correspondant a la fenetre 1
heat = G(H);                        % indices correspondant au radiateur
b(window1) = -1/h^2*ot;              % CL Dirichlet fenêtre 1
b(window2) = -1/h^2*ot;               % CL Dirichlet fenêtre 2
b(heat) = -ht;                      % CL Dirichlet chauffage
b(door) = -1/h^2*dt;                 % CL Dirichlet porte
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');                         % on dessine la solution sur la grille