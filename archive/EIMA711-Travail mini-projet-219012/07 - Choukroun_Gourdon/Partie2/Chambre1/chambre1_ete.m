% Ce script permet la simulation de refroidissement d'une pièce en été avec une
% température extérieure à 35°C, une température à la porte de 20°C, 
% et le chauffage (en mode climatiseur) apportant du froid à une puissance de 184 
% On utilise ici la méthode d'Euler explicite et implicite

ot = 35;
dt = 20;
ht = -184;
n = 30;

x = linspace(-1,1,n);   % x : vecteur de -1 à 1 avec n valeurs équidistantes
% X : matrice avec les valeurs de x en ligne donc toutes ses lignes sont égales à x
% Y : matrice avec les valeurs de x en colonnes donc toutes ses colonnes sont égales à x
[X,Y] = meshgrid(x,x);                                  % grille bi-dimensionnelle

% géométrie de la chambre
G = ((X>-1) & (X<1) & (Y>-0.3) & (Y<0.9)) | ...         % piece principale
    ((X>-0.7) & (X<-0.2) & (Y>-0.9) & (Y<-0.3)) | ...   % salle de bain
    ((X>0.5) & (X<0.8) & (Y>-0.9) & (Y<-0.3)) | ...     % dressing
    ((X>0.3) & (X<0.7) & (Y>0.9) & (Y<1)) | ...         % fenêtre du sud
    ((X>-0.55) & (X<-0.35) & (Y>-1) & (Y<-0.9));        % fenêtre du nord (salle de bain)
       
% position du chauffage
% bon positionnnement
% Chauffage 1 meilleur
H = ((X>0.2) & (X<0.7) & (Y>0.6) & (Y<0.9)) | ...        % chauffage en dessous de la fenêtre du mur sud
    ((X>-0.6) & (X<-0.3) & (Y>-0.9) & (Y<-0.65));        % chauffage salle de bain sur le mur nord

% mauvais positionnement
% Chauffage 2
%H = ((X>-0.8) & (X<-0.5) & (Y>0.7) & (Y<0.9)) | ...     % chauffage dans le coin ouest-sud
    %((X>-0.6) & (X<-0.3) & (Y>-0.9) & (Y<-0.6));        % chauffage salle de bain sur le mur nord

% Chauffage 3
%H = ((X>-0.25) & (X<0.25) & (Y>-0.2) & (Y<0.2))         % un grand chauffage au centre

% Chauffage 4
%H = ((X>0.35) & (X<0.65) & (Y>0.7) & (Y<0.9));          % chauffage uniquement en dessous de la fenêtre du mur sud

% Chauffage 5
%H = ((X>-0.8) & (X<-0.5) & (Y>0.7) & (Y<0.9));          % chauffage uniquement sur le mur ouest

k = find(G);                        % on trouve les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel % G garde sa taille mais tous ses éléments deviennent nuls 
G(k) = (1:length(k))';              % suite d'indices de la matrice A % G reprend sa forme de carré central avec cette fois les valeurs du carré allant de 1 au nombre de valeurs, s'implémentant dans le sens des colonnes 
A = delsq(G);                       % matrice du Laplacien à l'intérieur du domaine G
door = [];                          % indices correspondant à la porte
window1 = [];                       % indices correspondant à la fenetre du bas
window2 = [];                       % indices correspondant à la fenetre de la salle de bain (mur du haut)

for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              
        if (Y(i,j)>0.1 && Y(i,j)<0.4)     % Position de la porte (sur le mur ouest)
          door = [door ind];              % On garde la CL Dirichlet pour la porte 
        else
          A(ind,ind) = A(ind,ind)-1;      % CL Neumann sur le mur ouest
        end
      end
      if G(i,j+1) == 0                    % CL Neumann sur le mur est
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0                    % On garde la CL Dirichlet pour la fenêtre sud
          if (X(i,j)>0.3 && X(i,j)<0.7) 
          window1 = [window1 G(i,j)];     % indices correspodant à la fenêtre sud
          else
              A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur sud
          end
      end
      if G(i-1,j) == 0                    % On garde la CL Dirichlet pour la fenêtre nord
          if (X(i,j)>-0.55 && X(i,j)<-0.35) 
          window2 = [window2 G(i,j)];     % indices correspodant a la fenêtre nord
          else
              A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur nord
          end
      end 
    end
  end
end

h = 2/(n-1);
A = -A/h^2;                         % division par h^2
b = zeros(length(k),1);
heat = G(H);                        % indices correspondant au radiateur
b(window1) = -1/h^2*ot;             % prise en compte des CL Dirichlet
b(window2) = -1/h^2*ot;
b(heat) = -ht;                      % et du chauffage
b(door) = -1/h^2*dt;

[U,u] = Chambre1(ot,dt,0,n); % solution initiale
[U_final,u_final] = Chambre1(ot,dt,ht,n); % solution finale

% Paramètres de la discrétisation
L = 5 ;            % longueur du domaine
nx = 50 ;          % nombre de mailles
dx = L/(nx+1) ;    % pas d'espace
cfl = 0.1 ;        % cfl
nu = 1;            % coefficient de diffusion
dt2 = (dx^2/(2*nu))*cfl ; % dt = pas de temps % car cfl/2 = nu*dt/dx^2 car en dimension 2, c'est la somme des cfl pour chaque dimension
nt = 10 ;          % nombre de pas de temps effectues

[m,n]=size(A);
I = eye(m);

tic
while (mean(mean(u(G(G>0)))) > 22) % On boucle jusqu'à atteindre une température moyenne de la pièce de 22°C
    u_old = u;
    %u = (I+nu*A*dt2)*u_old-b*dt2; % explicite % On fait -b parce que dans l'équation de Poisson, on a +Laplacien de u = f et dans l'équation de la chaleur on on -nu*Laplacien de u = f
    u = (I-nu*dt2*A)\(u_old - nu*dt2*b); % implicite
    
    % U(G>0) et G(G>0) représentent les éléments supérieur strictement de U ou de G vu
    % qu'ils sont égaux - correspond à suite 1,2,3,...,length(k)
    % u(G(G>0)) correspond donc au k premiers éléments de u 
    % U(G>0) = full(u(G(G>0)) va faire de G une matrice gardant sa forme de
    % carré central mais au lieu d'une suite de 1 à length(k), on a les
    % solutions u de Au=b dans l'ordre en parcourant les colonnes
    U = G;
    U(G>0) = full(u(G(G>0)));  % on met la solution dans une matrice
    mesh(X,Y,U);               % correspondant à la grille
                               % De cette facon G a la taille de la grille et on peut donc la modéliser
                               % mesh va plot les valeurs de U dans le repère X,Y
    axis('ij');                % on dessine la solution sur la grille
    pause(0.001)               % pour afficher plus vite
end
toc
disp("t "+mean(mean(u(G(G>0)))));              % Affiche la température moyenne de la pièce (théorique)
disp("tfinal "+mean(mean(u_final(G(G>0)))));   % Affiche la température moyenne de la pièce (expérimentale)
