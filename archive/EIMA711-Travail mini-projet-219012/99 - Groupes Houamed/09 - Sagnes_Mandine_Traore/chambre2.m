function Chambre2(ot,dt,ht,n)
% Chambre2 calcule la température ambiante dans une chambre
% U = Chambre2(ot, dt, ht, n): prend la température extérieure ot,
% la température de la porte dt et la température du chauffage ht et la
% nombre de points de grille n, puis calcule la température ambiante de la
% chambre
x = linspace(-1,1,n);
y = linspace(-1,1,n);                        
[X,Y] = meshgrid(x,y);                             % grille bi-dimensionnelle 
G = ((X>-0.2) & (X<0.8) & (Y>-0.8) & (Y<0)) | ...  % geometrie de la chambre
   ((X>-0.8) & (X<0.8) & (Y>0) & (Y<0.8)) | ...
   ((X>-0.2) & (X<0.6) & (Y>0.8) & (Y<0.9));
H = ((X>-0.2) & (X<0.6) & (Y>0.75) & (Y<0.8));       % position du chauffage
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];
window=[];                          % indices correspondant a la porte et la fenetre
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0         
          A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur de gauche
      end
      if G(i,j+1) == 0                      
          A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur de droite
      end
      if G(i-1,j) == 0   
        if (X(i,j)>0.4 && X(i,j)<0.8) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte sur le mur du haut
        else
          A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du haut 
        end
      end
      if G(i+1,j) == 0      
        if (X(i,j)>-0.2 && X(i,j)<0.6)
          window = [window ind];    %% On garde la CL Dirichlet pour la porte sur le mur du bas
        else   
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du bas
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
b(door) = -1/h^2*dt;
b(heat) = -ht;                      % et du chauffage
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');                       % on dessine la solution sur la grille
end