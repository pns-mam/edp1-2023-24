function U = Fonction_RoomTemperatureP2(ot,dt,ht,n,C0,t)

% La fonction calcule la température ambiante dans la chambre 1 pour la partie 2 du projet 

% Température fenêtre : ot (Température extérieure)
% Température porte : dt 
% Température du chauffage : ht
% Nombre de points de la grille : n
% Condition Initial : CO (Température initiale de la chambre)
% Nombre d'itérations  pour la boucle en temps : t


% On dessine la chambre dans un premier temps avec les conditions de Neumann et de Dirichlet 
% puis on place le chauffage, la porte et la fenêtre
% on implémente une boucle pour le temps, pour voir l'évolution de la température dans la chambre
% et on résout grâce aux schémas implicite ou explicite l'équation de la chaleur 
% Et on affiche la température ambiante au cours du temps dans la chambre 1




x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             
G  = ((X>-1) & (X<-0.7) & (Y>-0.3) & (Y<0.3)) | ((X>-0.7) & (X<0.2) & (Y>-1) & (Y<1)) | ((X>0.2) & (X<1) & (Y>-1) & (Y<0)); % Dimension Chambre 1


H = ((X>-0.4) & (X<0.2) & (Y>-1) & (Y<-0.7)); % position du chauffage

k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A 

spy(G)                                % on dessine la grille
heat = G(H);                          % indices correspondant a H

 
A = delsq(G);                         % matrice du -Laplacien sans 1/h^2 a l'interieur
door = [];                            % indices correspondant a la porte
window = [];                          % indices correspondant a la fenetre
for i = 2:n-1                        
  for j = 2:n-1
    ind = G(i,j);
    if ind ~= 0 
      if G(i+1,j) == 0                % la porte se situe sur le mur du bas
        if (X(i,j)>-0.5 && X(i,j)<-0.1) 
          door = [door ind];          % indices correspondant a la porte 
        else
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la porte le mur est isolant
        end
      end
      if G(i-1,j) == 0                % mur du haut isolant: Neumann
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i,j+1) == 0                % mur a droite isolant: Neumann
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i,j-1) == 0                % la fenetre se situe sur le mur a gauche
        if (Y(i,j)>-0.2 && Y(i,j)<0.2) 
          window = [window G(i,j)];   % indices correspodant a la fenetre 
        else
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la fenetre le mur est isolant
        end
      end
    end
  end
end

h = 2/(n-1);                          % pas du maillage
A = -A/h^2;                           % division par h^2
b = zeros(length(k),1);               % initialisation du second membre
b(heat) = -ht;                        % prise en compte du chauffage
b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*ot;                % CL Dirichet: fenetre
alpha = 0.4;            % Condition CFL
 nu =1;                   % Coefficient de diffusion
deltat = alpha/(2*nu)*h^2;            %Pas de temps
u = C0*ones(length(k),1);               %Condition Initiale

tic   % Mesure le temps d'éxecution

for k = 1:t              % boucle en temps 
    uold = u;
    u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
    %u = (eye(size(A))+nu*deltat*A)*uold - nu*deltat*b;                            % explicite
    U = G;
    U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice
    % correspondant a la grille
    mesh(X,Y,U);                          % on dessine la solution sur la grille
    axis('ij');
    pause(0.001);
end
toc
