%% Partie 1

function U = RoomTemperature(ot,dt,ht,n)
% ROOMTEMPERATURE calcule la température ambiante dans une chambre
% U = RoomTemperature(ot, dt, ht, n); prend la température extérieure ot,
% la température de la porte dt et la température du chauffage ht et la
% nombre de points de grille n, puis calcule la température ambiante de la
% chambre

ot = 20; % temperature de la fenetre
dt = 20;  % temperature de la porte d'entree
ht = 0; % temperature du chauffageG
n = 30;   % nombre de points de discretisation


x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-0.7) & (X<1) & (Y>-0.4) & (Y<1)) |  ... % pièce principale
   ((X>0.7) & (X<1) & (Y>-1) & (Y<-0.4)) | ...    %couloir
   ((X>-1) & (X<-0.7) & (Y>-0.4) & (Y<0));        %petite extension
H = ((X>0.4) & (X<0.8) & (Y>0) & (Y<0.5)) | ...   % position du chauffage dans le couloir
    ((X>-0.7) & (X<-0.4) & (Y>-0.2) & (Y<0.2));    % position du chauffage dans la salle principale
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
      if G(i,j+1) == 0              % CL Neumann sur le mur droit
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0              % On garde la CL Dirichlet pour la fenetre
        if (X(i,j)>-0.2 && X(i,j)<0.3)
          window = [window ind];
        else
          A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas
        end 
      end
      if G(i-1,j) == 0
        if (X(i,j)>0.6 && X(i,j)<1) 
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
U(G>0) = full(u(G(G>0)));         % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');                         % on dessine la solution sur la grille

nb=0;
for i=1:length(k)
  if 18<= u(i) && u(i)<=23
    nb+=1;
  end
end
moyenne=nb/length(k)

max=max(u)

