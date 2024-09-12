ot=-15
dt=20
ht=0
n=40

% ROOMTEMPERATURE calcule la température ambiante dans une chambre
% U = RoomTemperature(ot, dt, ht, n); prend la température extérieure ot,
% la température de la porte dt et la température du chauffage ht et la
% nombre de points de grille n, puis calcule la température ambiante de la
% chambre

 x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-.9) & (X<-.4) & (Y>-0.8) & (Y<0)) | ...  % geometrie de la chambre
   ((X>-0.4) & (X<0.8) & (Y>-.8) & (Y<.4) & (Y-X<.4)) | ... 
   ((X>0) & (X<.8) & (Y>.4) & (Y<.9));
H = ((X>0.1) & (X<0.6) & (Y>0.6) & (Y<0.8));       % position du chauffage gauche
H = ((X>-0.7) & (X<0.5) & (Y>-0.6) & (Y<-.1));      % position du chauffage du bas
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = []; % indices correspondant a la porte
window = []; % indices correspondant a la fenêtre
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0  
        if (Y(i,j)>0) % On garde la CL Dirichlet pour la fenêtre de droite 
            window = [window ind]; 
        else % CL Neumann sur le mur gauche
            A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i,j+1) == 0      
        if (Y(i,j)>-0.8 && Y(i,j)<.5)  
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        else % CL Neumann sur le mur droit
          A(ind,ind) = A(ind,ind)-1; 
        end
      end
      if G(i+1,j) == 0 && i<n-1  
          if (X(i,j)<-.4)  % On garde la CL Dirichlet pour la fenêtre de gauche
            window = [window ind]; 
        else % CL Neumann sur le mur du bas
            A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i-1,j) == 0 % CL Neumann sur le mur du haut
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
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');                         % on dessine la solution sur la grille


