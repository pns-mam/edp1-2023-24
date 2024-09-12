clear all;
close all;
function U = RoomTemperature(ot,dt,ht,n)
% ROOMTEMPERATURE calcule la température ambiante dans une chambre
% U = RoomTemperature(ot, dt, ht, n); prend la température extérieure ot,
% la température de la porte dt et la température du chauffage ht et la
% nombre de points de grille n, puis calcule la température ambiante de la
% chambre

%Pour la partie stationnaire

%temperature ambiante en été ot=20,dt=20
%temperature ambiante en hiver ot=-10, dt=15
%Pour les deux premiere question de la question 2 on fixe ht à 0
%Pour la question ou le chauffage est allumé on fixe la puissance du chauffage à 110

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = (((X>-1) & (X<0.5) & (Y>0) & (Y<1))|...
     ((X>-0.2) & (X<0.5) & (Y>-1) & (Y<0.5))|...
     ((X>0.5) & (X<1) & (Y>-0.35) & (Y<0.2))) ;
     
%disp(G); si on veut afficher notre chambre sous forme de matrice 
%spy(G); si on veut afficher notre chambre dans un repère comme dans notre rapport


H = ((X>0) & (X<0.5) & (Y>0.5) & (Y<1));      % chauffage au centre

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
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i,j+1) == 0 && j<n-1         % CL Neumann sur le mur droit
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0      % On garde la CL Dirichlet pour la fenetre
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
      end
      if G(i-1,j) == 0
        if (X(i,j)>-1 && X(i,j)<-0.2) 
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
window = G(G(:,end-1)>0,end-1);     % indices correspondant a la fenetre

heat = G(H);                        % indices correspondant au radiateur faire attention qu'il se trouve bien dans la pièce

b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(heat) = -ht;                      % et du chauffage
b(door) = -1/h^2*dt;
I = eye(size(A));
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij')
end

U = RoomTemperature(-10,15,138,40);  % Permet de définir les valeurs initiales de notre fonction
x = linspace(-1,1,40);                            
[X,Y] = meshgrid(x,x);
MoyU = mean(mean(U(((X>-1) & (X<0.5) & (Y>0) & (Y<1))|...
     ((X>-0.2) & (X<0.5) & (Y>-1) & (Y<0.5))|...
     ((X>0.5) & (X<1) & (Y>-0.35) & (Y<0.2)))));     

disp(MoyU); % On affiche la valeur moyenne de notre chambre dans le terminal
