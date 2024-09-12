% Partie 2 : Problème instationnaire

function U = RoomTemperature_instationnaire_PIERROT(wt,dt,ht,n)
% Probleme de Dirichlet non-homogene
% Avec deux termes source localises (de type chauffage ou climatisation)
% Prise en compte des CL Dirichlet (ajout d'une fenetre + porte)
wt = -10; % temperature de la fenetre
dt = 10;  % temperature de la porte d'entree
ht = 300; % temperature du chauffage
n = 30;   % nombre de points de discretisation

x = linspace(-1,1,n);
[X,Y] = meshgrid(x,x);                  % grille bi-dimensionnelle
G = ((X>-0.7) & (X<1) & (Y>-0.4) & (Y<1)) |  ... % pièce principale
   ((X>0.6) & (X<1) & (Y>-1) & (Y<-0.4)) | ...    %couloir
   ((X>-1) & (X<-0.7) & (Y>-0.4) & (Y<0));        %petite extension sur la gauche
H = ((X>0.4) & (X<0.8) & (Y>0) & (Y<0.5)) | ...   % position du chauffage dans le couloir
    ((X>-0.7) & (X<-0.4) & (Y>-0.2) & (Y<0.2));    % position du chauffage dans la salle principale


k = find(G);                          % les indices des G non nuls
G = zeros(size(G));                   % conversion logique - reel
G(k) = (1:length(k))';                % suite d'indices de la matrice A
spy(G)                                % on dessine la grille
heat = G(H);                          % indices correspondant a H


A = delsq(G);                         % matrice du Laplacien a l'interieur
door = [];                            % indices correspondants à la porte
window = [];                          % indices correspondant à la fenetre


for i = 2:n-1
    for j = 2:n-1
        ind=G(i,j);
        if G(i,j) ~= 0
            if G(i+1,j) == 0                % la fenetre se situe sur le mur du bas
                if ((X(i,j)>-0.2 && X(i,j)<0.3))
                    window = [window ind];   % indices correspodant a la fenetre
                else
                    A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas
                end
            end
            if G(i-1,j) == 0                % la porte se situe sur le mur du haut
              if (X(i,j)> 0.6 && X(i,j)< 1)
                  door = [door ind];       % indices correspondant a la porte
              else
                  A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du haut
              end  
            end
            if G(i,j+1) == 0                
                A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur de droite
            end
            if G(i,j-1) == 0                
                A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur de gauche
            end
        end
    end
end


h = 2/(n-1);                          % pas du maillage
A = -A/h^2;                           % division par h^2
b = zeros(length(k),1);               % initialisation du second membre
b(heat) = -ht;                        % prise en compte du chauffage
b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*wt;                % CL Dirichet: fenetre


cfl = 0.4, nu =1;
deltat = cfl/(nu)*h^2;
u = 30*ones(length(k),1);     % solution initiale homogène en température
iteration=100                 % nombre d'itération minimum permettant une quasi-convergence
tic                            % début du chronomètre
for l = 1:iteration                % boucle en temps 
    uold = u;
    %u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
    %u = (eye(size(A))+nu*deltat*A)*uold - nu*deltat*b; % explicite
    u=(eye(size(A))-nu*deltat*A)\((eye(size(A))+nu*deltat*A)*uold-nu*deltat*b); %crank-nicholson
    u(door) = dt;                                     %temperature de la porte fixée                
    u(window) = wt;                                   % température fenetre fixée
    U = G;
    U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice
    % correspondant a la grille
    %mesh(X,Y,U);                          % on dessine la solution sur la grille
    %axis('ij');
    %pause(0.05)
end
toc                                   % fin du chronomètre
mesh(X,Y,U);                          % on dessine la solution sur la grille
axis('ij');


nb=0;
for i=1:length(k)
  if 18<= u(i) && u(i)<=23    % calcul du pourcentage de la chambre considéré 
    nb+=1;                    % etant confortable en terme de temperature              
  end                         % (entre 18°C et 23°C)
end
moyenne=nb/length(k)           

max=max(u)                    % calcul de la température maximale
min=min(u)                    % calcul de la température minimale
