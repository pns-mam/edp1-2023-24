
tic
ot=-15;     % temperature exterieure donc celle de le fenetre
dp=15;     % la temperature de la porte dp
ht=600;   %la temperature du chauffage ht
n=40;      % nombre de points de grille n

% Nous allons calculer la température ambiante dans la chambre 1
% a l'aide de l'equation de la chaleur en mode instationnaire.
% Nous utilisons le schéma explicite ou le schéma implicite pour nos
% simulations. Ici nous simulons l'évolution dans la chambre 1 froide
% à partir du moment où on allume le chauffage.


x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
   
G = ((X>-1) & (X<0.6) & (Y>0) & (Y<0.8)) | ...  % geometrie de la chambre
   ((X>-1) & (X<0.6) & (Y>-0.5) & (Y<0)) | ...
   ((X>-0.4) & (X<0.4) & (Y>-1) & (Y<-0.5)) | ...
   ((X>0.6) & (X<0.75) & (Y>-0.3) & (Y<0.3));

H = ((X>0.3) & (X<0.55) & (Y>-0.3) & (Y<0.3));       % position du chauffage
                                                     % en dessous de la fenetre
%H = ((X>-0.3) & (X<0.3) & (Y>0.5) & (Y<0.75));      % position du chauffage
                                                     % sur le mur au fond
%H = ((X>-0.95) & (X<-0.7) & (Y>-0.3) & (Y<0.3));       % position du chauffage sur le mur de gauche
                                                     
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
window = [];                        % indices correspondant a la fenÃªtre
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              % CL Neumann sur le mur gauche
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0      
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
      end
      if G(i,j+1) == 0                % la fenetre se situe sur le mur a droite
        if (Y(i,j)>-0.3 && Y(i,j)<0.3) 
          window = [window G(i,j)];   % indices correspodant a la fenetre 
        else
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la fenetre le reste du mur est isolant
        end
      end
      if G(i-1,j) == 0
        if (X(i,j)>-1 && X(i,j)<-0.5) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        else
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
        end
      end
    end
  end
end

L = 2;            % longueur du domaine
dx = L/(n-1) ;    % pas d'espace
cfl = 0.5 ;        % condition cfl
nu = 1;            % coefficient de diffusion
dt = cfl/(2*nu)*dx^2; % dt = pas de temps optimal pour le schéma explicite
nt = 3000 ;          % nombre d'itérations maximales

A = -A/dx^2;       % Calcul du Laplacien à l'interieur de la grille G

b = zeros(length(k),1);
heat = G(H);                        % indices correspondant au radiateur
b(window) = -1/dx^2*ot;             % prise en compte des CL Dirichlet pour la fenêtre
b(heat) = -ht;                      % prise en compte du chauffage 
b(door) = -1/dx^2*dp;               % prise en compte des CL Dirichlet pour la porte 
u0 = 10*ones(length(k),1);          % initialisation de la solution
u0(door) = dp;                      % prise en compte du de l'energie du chauffage dans la solution initiale 
u0(window) = ot;                    % prise en compte du froid de la fenetre la solution initiale
temp = sum(u0)/length(k);    % calcul de la température ambiante a t=0 en faisant une moyenne
it = 0;                      % varible pour itérer

% On utilise une boucle while pour calculer la solution jusqu'a ce que
% la température demandée soit atteint ou qu'on dépasse le nbr max
% d'itérations
while (temp<20 && it<nt)
    %u = (eye(length(u0),length(u0))-nu*dt*A)\(u0 - nu*dt*b); % calcul avec schema implicite
    u = (eye(length(u0),length(u0))+nu*dt*A)*u0-nu*dt*b; % calcul avec shema explicite
    
    u0 = u;  % on stocke la solution
    pause(0.01)
    
    U = G;
    U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
    mesh(X,Y,U);                        % correspondant a la grille 
    axis('ij');                         % on dessine la solution sur la grille
    title('Chauffage au niveau de la fenetre');
    temp = sum(sum(U))/length(k); % calcul de la température ambiante au temps actuel
    it = it+1; % on incremente
end
toc


