clear all
% Simulation du probleme de Poisson
% Probleme de Dirichlet non-homogene
% Avec un terme source localise (de type chauffage)
% Prise en compte des CL Dirichlet (ajout d'une fenetre + porte)

function instationnaire(n,ot,ht,dt,u0)
tic                                                %debut temps d execution
x = linspace(-1,1,n);
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle

%Géométrie des elements a ajouter
G = ((X>-1) & (X<-0.4) & (Y>-1) & (Y<0.8)) | ...       %geometrie de la piece
    ((X>-0.4) & (X<0.7) & (Y>0.3) & (Y<0.8))|  ... 
    ((X>-0.1) & (X<0.3) & (Y>0.8) & (Y<1)); 
H = ((X>-0.1) & (X<0.4) & (Y>0.5) & (Y<0.7));         %chauffage fenetre

k = find(G);                          % les indices des G non nuls
G = zeros(size(G));                   % conversion logique - reel
G(k) = (1:length(k))';                % suite d'indices de la matrice A
spy(G);                               % on dessine la grille
heat = G(H);                          % indices correspondant a H

A = delsq(G);                         % matrice du Laplacien a l'interieur

%boucles d ajout des conditions de Neumann pour les murs isolants, portes et fenetres
door = [];
window = [];
for i = 2:n-1                       
  for j = 2:n-1                     
    ind = G(i,j);
    if ind~=0
      if G(i,j+1) == 0              % CL Neumann sur le mur droit
         A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0
        if ((X(i,j)>-0.1) && (X(i,j)<0.3)) % indices correspondant a la porte
          window = [window ind];    % On garde la CL Dirichlet pour la fenetre
        else
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
        end
      end
      if G(i-1,j) == 0
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
      end
      if G(i,j-1) == 0              % CL Neumann sur le mur gauche
          if (Y(i,j)>-0.5 && Y(i,j)<0) 
            door = [door ind];        % On garde la CL Dirichlet pour la porte 
          else
            A(ind,ind) = A(ind,ind)-1;
          end
      end
    end
  end
end


h = 2/(n-1);                          % pas du maillage
A = -A/h^2;                           % division par h^2
b = zeros(length(k),1);               % initialisation du second membre


%ajout des conditions limites pour la porte, fenetre et chauffage
b(heat) = -ht;                        % prise en compte du chauffage
b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*ot;                % CL Dirichet: fenetre


%parametre de la discrétisation
nu =1;                                %coefficient de diffusion
L=5;                                  %longueur du domaine
dx = L/(n+1);                         %pas d'espace
cfl = 0.04;                           %cfl
delta = dx^2/nu*cfl;                  %pas de temps

%ajout de notre boucle en temps
u = u0*ones(size(b));                 %temperature dans la salle a l'initialisation
v = u0*ones(size(b));                 %temperature dans la salle a l'initialisation

for l = 1:30                          %boucle en temps 
    uold = u;
    %u = (eye(size(A))-nu*delta*A)\(uold - nu*delta*b); %euler implicite
    u = (eye(size(A))+nu*delta*A)*uold - nu*delta*b;      %euler explicite
    U = G;
    U(G>0) = full(u(G(G>0)));          %on met la solution dans une matrice
               % correspondant a la grille

    mesh(X,Y,U);                       %on dessine la solution sur la grille
    axis('ij');
    title(['t= ', num2str(l)]);
    pause(0.1);
end

%res=V-U;                              %comparaison entre explicite et implicite
%disp(res);

toc                                    %fin du temps d execution
endfunction


%realisation de nos cas tests instationnaire pour euler explicite
%pour utiliser euler implicte mettre en commentaire euler explicite et decommenter euler implicite
%n maillage
%ot temperature exterieur
%ht temperature chauffage
%dt temperature porte
%uO temperature initiale dans la piece

%cas 1 : très chaud dehors et ajout de la climatisation
instationnaire(30,30,-200,20,30);

%cas 2 : très froid dehors et ajout du chauffage
%instationnaire(30,-10,400,5,10);




