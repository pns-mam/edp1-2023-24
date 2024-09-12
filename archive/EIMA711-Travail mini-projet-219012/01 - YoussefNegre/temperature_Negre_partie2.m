clear all, close all, clc



% ---- VARIABLES ----
  ht = 70
  dt = 25
  wt = 0
  n = 50
  iteration = 20
  ini = 10
  methode = 0
  test = 1
% Simulation du l'equation de la chaleur.
% Probleme mixte: de Dirichlet non-homogene + Neumann homogene 
% Terme source localise (de type chauffage)
% Prise en compte des CL Dirichlet (ajout de 3 fenetres + porte)
% Prise en compte des CL Neumann (murs isolants)
%ht  temperature du chauffage
%dt  temperature de la porte d'entree
%wt  temperature de la fenetre
%n   nombre de points de discretisation
%iteration nombre d'iteration dans la boucle temporelle
%ini température initiale dans l'ensemble de la pièce (temperature homogène) 
%methode type de methode utilisée : explicite = 0 / implicite = 1  
%test fonction du radiateur : 0 = eteint, 1 = chauffage progressif, 2 = climatiseur , 3 = chauffage constant

x = linspace(-1,1,n);                      
[X,Y] = meshgrid(x,x);                  % grille bi-dimensionnelle 

G = ((X>-1) & (X<-0.2) & (Y>-1) & (Y<1)) | ...  % geometrie de la chambre
   ((X>-0.2) & (X<0.5) & (Y>0.6) & (Y<1)) | ... 
   ((X>-0.2) & (X<0.7) & (Y>-0.5) & (Y<-0.2)) | ((X>-0.2) & (X<0.2) & (Y>-1) & (Y<-0.7)) ;
   
H = ((X>-1) & (X<-0.6) & (Y>-0.9) & (Y<-0.6));       % position du chauffage

k = find(G);                          % les indices des G non nuls
G_ = G;
G = zeros(size(G));                   % conversion logique - reel
G(k) = (1:length(k))';                % suite d'indices de la matrice A
spy(G)                                % on dessine la grille
heat = G(H);                          % indices correspondant a H
A = delsq(G);                         % matrice du -Laplacien sans 1/h^2 a l'interieur;
door = [];                            % indices correspondant a la porte
window = [];                          % indices correspondant a la fenetre
for i = 2:n-1                        
  for j = 2:n-1
    ind = G(i,j);
    if ind ~= 0 
      if G(i-1,j) == 0                % la porte se situe sur le mur du haut
        if (X(i,j)>-0.4 && X(i,j)<-0.3) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
         else
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la porte le mur est isolant
        end
      end
      if G(i+1,j) == 0                % mur du bas isolant: Neumann
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i,j-1) == 0                % mur a gauche isolant: Neumann
        if (Y(i,j)>-0.4 && Y(i,j)<0.3) 
          window = [window G(i,j)];   % indices correspodant a la fenetre 
        else
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i,j+1) == 0                % la fenetre se situe sur le mur a droite
        
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la fenetre le mur est isolant
        end
      end
    end
  end
end

h = 2/(n-1);                          % pas du maillage
A = -A/h^2;                           % division par h^2
b = zeros(length(k),1);               % initialisation du second membre

if test != 0
  b(heat) = -ht;                      % prise en compte du chauffage  
end

b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*wt;                % CL Dirichet: fenetre

alpha = 0.5; nu =1;
deltat = alpha/(2*nu)*h^2;
u = ini*ones(length(k),1) + zeros(length(k),1);     % solution initiale 




tic()
for l = 1:iteration         % boucle en temps 
    uold = u;
    
    if methode == 0
          u = (eye(size(A))+nu*deltat*A)*uold - nu*deltat*b; % explicite
    end
    if methode == 1
        u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
    end

    
    if test == 1              %Le radiateur augmente en temperature
      if u(heat)< ht-0.5 
       u(heat) = u(heat) +0.5;
      else
       u(heat) = ht;
      end
    end
    if test == 2              %Le climatiseur diminue en temperature
      if u(heat)> ht +0.5
       u(heat) = u(heat) - 0.5;
      else
       u(heat) = ht;
      end
    end
    if test == 3            %Le radiateur de temperature ht à tout instant.
      u(heat) = ht;
    end
    U = G;
    U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice
    
    
    % correspondant a la grille
     mesh(X,Y,U);                          % on dessine la solution sur la grille
    axis([-1 1 -1 1 0 80],"ij");
    title({"t = ";l*deltat;"Itération : ";l});
    pause(0.001)
end 
 toc() 

 
