%{
MINI PROJET - EDP1 MAM4
Mod�lisation de la temp�rature dans une pi�ce
DOUGNAC Jade / BARNETCHE Carine
%}

% PARTIE 2 : CHAMBRE 1 - MODELE INSTATIONNAIRE

clear all, close all, clc
function U = RoomTemperature(ot,dt,ht,n)
%{ 
ROOMTEMPERATURE calcule la temp�rature ambiante dans une chambre
U = RoomTemperature(ot, dt, ht, n); prend la temp�rature ext�rieure ot,
la temp�rature de la porte dt et la temp�rature du chauffage ht et la
nombre de points de grille n, puis calcule la temp�rature ambiante de la
chambre 
%}
  x = linspace(-1,1,n);                            
  [X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
  G = ((X>0.4) & (X<0.8) & (Y>-1) & (Y<-0.5)) | ...  % geometrie de la chambre
      ((X>-1) & (X<-0.3) & (Y>-0.6) & (Y<0.2)) | ... 
      ((X>-0.3) & (X<0.8) & (Y>-0.5) & (Y<0.5))|...
      ((X>-0.7) & (X<0.2) & (Y>0) & (Y<0.8));

  H = ((X>-0.6) & (X<0) & (Y>-0.3) & (Y<0));          % chauffage 

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
        if G(i,j-1) == 0              
          A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur gauche
        end
        if G(i,j+1) == 0            
          A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur droit
        end
        if G(i+1,j) == 0 
          if (X(i,j)>-0.4 && X(i,j)<0)
            window = [window ind];    % On garde la CL Dirichlet pour la fenetre
          else
            A(ind,ind) = A(ind,ind)-1;% CL Neumann sur le mur du bas 
          end
        end
        if G(i-1,j) == 0 
          if (X(i,j)>-0.7 && X(i,j)<-0.4) 
            door = [door ind];        % On garde la CL Dirichlet pour la porte 
         else
            A(ind,ind) = A(ind,ind)-1;% CL Neumann sur le mur du haut
          end
        end 
      end
    end
  end


  h = 2/(n-1);                        % pas d'espace
  A = -A/h^2;                         % division par h^2
  b = zeros(length(k),1);

  heat = G(H);                        % indices correspondant au radiateur

  b = zeros(length(k),1);             % initialisation du second membre
  b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
  b(door) = -1/h^2*dt;
  b(heat) = -ht;                      % et du chauffage

  alpha = 0.5; nu = 1;
  deltat = alpha/(2*nu)*(h^2);
  I = eye(size(A));
  u = ot*ones(length(k),1);  % au d�part la pi�ce est � la temp�rature ext�rieure

  tic 
  for i = 1:3848                        % boucle en temps
      uold = u;
      %u = (I - nu*deltat*A)\(uold - nu*deltat*b); % implicite
      u = (I + nu*deltat*A)*uold - nu*deltat*b;    % explicite
      U = G;
      U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
      %{ 
      Pour obtenir les temps d'executions on met en commentaire 
      la boucle if qui est ici pour l'affichage 
      %}
      if mod(i,100) == 0
        mesh(X,Y,U);                      % correspondant a la grille 
        axis('ij');                       % on dessine la solution sur la grille
        pause(0.1);                     
      end     
  end
  toc
  
end

U = RoomTemperature(-10,-10,500,40); 

% On affiche la temp�rature finale moyenne dans la pi�ce 
% on fait bien attention � calculer les valeurs � l'int�rieur de la chambre
x = linspace(-1,1,40);                            
[X,Y] = meshgrid(x,x); 
temperature_moyenne = mean(mean(U((((X>0.4) & (X<0.8) & (Y>-1) & (Y<-0.5)) | ...  
    ((X>-1) & (X<-0.3) & (Y>-0.6) & (Y<0.2)) | ... 
    ((X>-0.3) & (X<0.8) & (Y>-0.5) & (Y<0.5))|...
    ((X>-0.7) & (X<0.2) & (Y>0) & (Y<0.8))))));
disp(temperature_moyenne);
