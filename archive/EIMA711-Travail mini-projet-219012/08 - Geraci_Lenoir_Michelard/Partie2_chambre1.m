% ------------- Cas Tests ------------- %
% Cas Test 1 : Hiver avec chauffage :
ot = -10;
dt = 15; 
ht = 600;
n = 40 ;
%% Cas Test 2 : Ete avec clim :
%ot = 30;
%dt = 20; 
%ht = -300;
%n = 40 ;


% ------------- Definition de l'espace (chambre 1) ------------- %
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-1) & (X<0.5) & (Y>-0.8) & (Y<0)) | ... 
    ((X>-0.7) & (X<0) & (Y>-1) & (Y<-0.8)) | ...
    ((X>-1) & (X<-0.25) & (Y>0) & (Y<0.5)) | ...
    ((X>-1) & (X<0) & (Y>0.5) & (Y<1));
   
%H = ((X>-0.5) & (X<0) & (Y>-0.8) & (Y<-0.3));     % chauffage au centre
H = ((X>-1) & (X<-0.75) & (Y>-0.8) & (Y<-0.3));      % chauffage près le mur
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = []; % indices correspondant a la porte
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
      if G(i+1,j) == 0 
        if (X(i,j)>-0.8 && X(i,j)<-0.3) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        else
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
        end
      end
      if G(i-1,j) == 0  % On garde la CL Dirichlet pour la fenetre
        if (X(i,j)>-0.7 && X(i,j)<0 && Y(i,j)>-1 && Y(i,j)<-0.8)
            window = [window ind];
        else
            A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
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


I = eye(size(A)); 

% solution initiale pour l'hiver :
u = ones(length(k),1)*5;  
u(door) = dt ; 
u(window) = 0 ; 

% solution initiale pour l'été :
%u = ones(length(k),1)*20;
%u(door) = dt ; 
%u(window) = ot ; 

nt = 1500;
alpha = 0.5; %cfl
nu =1;
deltat = alpha/(2*nu)*h^2;
t=0;
% On boucle sur le temps
tic
for t = 1:nt       
    uold = u;
    u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
    %u = (eye(size(A)) + nu * deltat * A) * uold - nu * deltat * b;% explicite
    U = G;
    U(G>0) = full(u(G(G>0)));% on met la solution dans une matrice correspondant a la grille
    % correspondant a la grille
    mesh(X,Y,U);
    title(t)
    axis('ij');
    pause(0.05);
end
toc

moy = sum(sum(U))/length(k); % On somme les colonnes puis le resultat des colonnes
disp("La temperature moyenne dans la chambre est de : ");disp(moy);