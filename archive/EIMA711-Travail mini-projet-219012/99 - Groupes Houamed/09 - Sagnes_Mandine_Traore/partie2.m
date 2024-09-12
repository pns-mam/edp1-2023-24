%{
  Ce fichier correspond à la réponse de la partie 2.
  Simulation stationnaire de l'équation de la chaleur, en deux dimensions.

  Définition des températures.
  Ce sont ces paramètres avec lesquels nous allons jouer pour effectuer différents tests...
  Définition du paramètre de discrétisation spatiale
%}

ht = 0; % temperature du chauffage
dt = 15;  % temperature de la porte d'entree
wt = -10; % temperature de la fenetre
n = 40;   % nombre de points de discretisation de notre grille
% flag qui correspond au choix de la géométrie de notre pièce
piece=2; % 1 (resp. 2) pour la géométrie de la chambre 1 (voir fichier Chambre1.m) (resp. chambre 2 (voir fichier Chambre2.m))
% flag qui correspond au choix de la méthode
methode=2; %1 (resp. 2) pour Euler explicite (resp. Euler implicite)

%{
  Dans cette section, nous allons définir
    la géométrie de la pièce (nous avons repris celles de la partie 1),
    stocker les indices du chauffage, porte, fenêtre
    prendre en compte Cl de Neumann en modifiant le laplacien là où il faut
      (si on est sur un bord du domaine on soustrait -1 au terme adéquat sur la diagonale de la matrice A)
%}
 
x = linspace(-1,1,n);
[X,Y] = meshgrid(x,x);                  % grille bi-dimensionnelle

if piece==1 
  G = ((X>0.2) & (X<0.8) & (Y>-0.6) & (Y<0.2)) | ...  % geometrie de la chambre 1
   ((X>-0.6) & (X<0.2) & (Y>-0.2) & (Y<0.8)) | ... 
   ((X>-0.7) & (X<-0.6) & (Y>0.2) & (Y<0.6));
  H = ((X>-0.6) & (X<-0.55) & (Y>0.2) & (Y<0.6));     % position du chauffage
end

if piece==2
  G = ((X>-0.2) & (X<0.8) & (Y>-0.8) & (Y<0)) | ...  % geometrie de la chambre 2
   ((X>-0.8) & (X<0.8) & (Y>0) & (Y<0.8)) | ...
   ((X>-0.2) & (X<0.6) & (Y>0.8) & (Y<0.9));
  H = ((X>-0.2) & (X<0.6) & (Y>0.75) & (Y<0.8));       % position du chauffage
end

k = find(G);                          % les indices des éléments de G non nuls (les points où il faudra calculer une température)
G = zeros(size(G));                   
G(k) = (1:length(k))';                % les éléments non-nuls de G (donc tous les 1) sont indexés(remplacés) par des entiers
A = delsq(G);                         % matrice du Laplacien a l'interieur (ici, on obtient l'opposé du laplacien théorique)
door = [];
window = [];

if piece==1
  for i = 2:n-1   %on parcourt les lignes
    for j = 2:n-1   %on parcourt les colonnes
      ind = G(i,j);   
      if ind~=0     %si on est dans la pièce
        if G(i,j-1) == 0       %si la colonne précédente est nulle      
          if (Y(i,j)>0.2 && Y(i,j)<0.6)  
            window = [window G(i,j)];    % On garde les indices où la CL Dirichlet doit être prise en compte, pour la fenêtre qui est sur le mur gauche
          else
            A(ind,ind) = A(ind,ind)-1;  % CL Neumann: pour un mur gauche isolant
          end
        end
        if G(i,j+1) == 0       %si la colonne suivante est nulle         
          A(ind,ind) = A(ind,ind)-1;    % CL Neumann: pour un mur droit isolant
        end
        if G(i-1,j) == 0       %si la ligne précédente est nulle  
          if (X(i,j)>0.2 && X(i,j)<0.6) 
            door = [door ind];         % On garde les indices où la CL Dirichlet doit être prise en compte, pour la porte qui est sur le mur du haut
          else
            A(ind,ind) = A(ind,ind)-1;  % CL Neumann: pour un mur du haut isolant 
          end
        end
        if G(i+1,j) == 0        %si la ligne suivante est nulle  
            A(ind,ind) = A(ind,ind)-1; % CL Neumann: pour un mur du bas isolant
        end 
      end
    end
  end
end

if piece==2
  for i = 2:n-1       % ajout de CL Neumann pour les murs isolants
    for j = 2:n-1                      
      ind = G(i,j);
      if ind~=0
        if G(i,j-1) == 0         
            A(ind,ind) = A(ind,ind)-1;  % CL Neumann: pour un mur gauche isolant
        end
        if G(i,j+1) == 0                      
            A(ind,ind) = A(ind,ind)-1;  % CL Neumann: pour un mur droit isolant
        end
        if G(i-1,j) == 0   
          if (X(i,j)>0.4 && X(i,j)<0.8) 
            door = [door ind];        % On garde les indices où la CL Dirichlet doit être prise en compte, pour la porte qui est sur le mur du haut
          else
            A(ind,ind) = A(ind,ind)-1;  % CL Neumann: pour un mur du haut isolant 
          end
        end
        if G(i+1,j) == 0      
          if (X(i,j)>-0.2 && X(i,j)<0.6)
            window = [window ind];    % On garde les indices où la CL Dirichlet doit être prise en compte, pour la fenetre qui est sur le mur du haut
          else   
            A(ind,ind) = A(ind,ind)-1; % CL Neumann: pour un mur du bas isolant
          end
        end 
      end
    end
  end
end

h = 2/(n-1);                          % pas spatial du maillage (deltaX = deltaY) selon les deux axes X et Y
A = -A/h^2;                           % On se ramène au laplacien théorique en divisant par -h^2 (à ce stade la matrice contient les CL de Neumann relativement à la matrice G)
heat = G(H);                          % on récupère les indices correspondant au chauffage


b = zeros(length(k),1);               % initialisation du membre qui contiendra les Cl de dirichlets aux indices adéquats
b(heat) = -ht;                        % prise en compte du chauffage
b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*wt;                % CL Dirichet: fenetre

%{
  Dans cette section nous allons, implémenter le schéma de résolution:
    déclarer les paramètres physiques (CFL, deltaT ...)
    déclarer une condition initiale ( que nous allons faire varier selon les tests...)
    implémenter la boucle temporelle 
%}

alpha = 1000000; % valeur du CFL
nbIteration=1000; %nombre de pas temporelle 

if (methode==1 && alpha>1/2) % si Euler explicite est choisi et que le CFL est supérieur à 0.5
  error("attention la schéma de résolution va être instable")
end

nu =1; %coefficient de diffusion
deltat = alpha/(2*nu)*h^2; %notre pas temporelle du schéma de résolution

%Condition initiale
%u = A\b;     % on part de la solution de l'équation de poisson
u = zeros(length(k),1)-10; % chambre initialement froide
%u = zeros(length(k),1)+25; % chambre initialement chaude

%on dessine la condition initiale
U = G;
U(G>0) = full(u(G(G>0)));             % on convertit la solution (les températures) en une matrice correspondant aux indices de la pièce
mesh(X,Y,U);                          % on dessine la solution sur la grille
title(strcat("At time",num2str(0)));
axis('ij');
pause(0.1)

tic
for l = 1:nbIteration             % boucle temporelle 
    uold = u;
    if methode==1 % explicite
      u = (eye(size(A))+nu*deltat*A)*uold - nu*deltat*b; 
    end
    if methode==2 % implicite
      u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); 
    end
    
    %on dessine la solution u
    U = G;
    U(G>0) = full(u(G(G>0)));             % on convertit la solution (les températures) en une matrice correspondant aux indices de la pièce
    mesh(X,Y,U);                          % on dessine la solution sur la grille
    title(strcat("At time ",num2str(l)));
    axis('ij');
    %pause(0.1)
    
%    %Si on veut augmenter la clime
%    if b(heat)(1)<=1500
%      b(heat)=b(heat)+100;
%    end
%    % on calcule la temperature ambiante (moyenne de la temp dans la pièce)
%    temp_ambiante=sum(u)/size(u)(1); 
%    if temp_ambiante<7.5
%      disp("on sort")
%      break
%    end
  
    %Si on veut augmenter le chauffage
    if b(heat)(1)>=-1500
      b(heat)=b(heat)-100;
    end
    % on calcule la temperature ambiante (moyenne de la temp dans la pièce)
    temp_ambiante=sum(u)/size(u)(1);
    if temp_ambiante>7.5
      disp("on sort")
      break
    end
end
toc