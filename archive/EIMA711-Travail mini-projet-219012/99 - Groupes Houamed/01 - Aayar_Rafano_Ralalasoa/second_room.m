%AAYAR - RAFANO - RALALASOA AVEC L'AIDE DE MME.DOLEAN
%%SCRIPT qui calcule la temp�rature ambiante dans une chambre
%%PARTIE 1 ET 2 INCLUSES

n=30;         %nombre de points de la grille
wt = -10;      %temp�rature ext�rieure (celle pr�s des fen�tres)
dt = 15;       % temp�rature de la porte
ht = 500;       % temp�rature du chauffage
  
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-1) & (X<-0.5) & (Y>-0.5) & (Y<0.5)) | ...  % geometrie de la chambre
   ((X>-0.5) & (X<0.5) & (Y>-1) & (Y<1)) | ... 
   ((X>0.5) & (X<1) & (Y>-0.5) & (Y<0.5));  
H = ((X>-0.2) & (X<0.2) & (Y>0.9) & (Y<1)) | ... % position des chauffages
    ((X>-0.2) & (X<0.2) & (Y>-1) & (Y<-0.9))
%P = ((X>0.9) & (X<1) & (Y>-0.2) & (Y<0.2))  %position de la porte 
%F = ((X>-1) & (X<-0.9) & (Y>-0.3) & (Y<0.3))   %position de la fen�tre
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  

%spy(G,"b");hold on;spy(H,"r");hold on; spy(P,"y");hold on;spy(F,"g") %pour v�rifier la position 


A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
window = [];                        %indices correspondant a la fen�tre
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0 && j<n-1             % CL Neumann sur le mur gauche
         if (Y(i,j)>-0.3 && Y(i,j)<0.3)  % On garde la CL Dirichlet pour la fen�tre
          window = [window ind];
          else
        A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i,j+1) == 0
        if (Y(i,j)>-0.2) && (Y(i,j)<0.2) % On garde la CL Dirichlet pour la porte
          door = [door ind];
          else
        A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur droit
        end
      end
      if G(i+1,j) == 0     
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas
      end
      if G(i-1,j) == 0 
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
      end 
    end
  end
end

%%%%%%%%%%%%%%%%PARTIE 1%%%%%%%%%%%%%%%%

h = 2/(n-1);
A = -A/h^2;                         % division par h^2
b = zeros(length(k),1);      % indices correspondant a la fenetre
heat = G(H);                        % indices correspondant au radiateur
b(heat) = -ht;                      % prise en compte du chauffage
b(window) = -1/h^2*wt;              % CL Dirichlet : fenetre
b(door) = -1/h^2*dt;                % CL Dirichlet : porte
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');                         % on dessine la solution sur la grille




%%%%%%%%%%%%%%%PARTIE 2%%%%%%%%%%%%%%%%%%%
%{
h = 2/(n-1);                          % pas d'espace
A = -A/h^2;                           % division par h^2

alpha = 0.5; nu =1;                   % CFL
deltat = alpha/(2*nu)*h^2;
heat = G(H);                          % indices correspondant au radiateur
b = zeros(length(k),1);               % initialisation du second membre
b(heat) = -ht;                        % prise en compte du chauffage
b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*wt;                % CL Dirichet: fenetre
uold = zeros(length(k),1);            % solution initiale pour cas test ski
%uold = 35*ones(length(k),1);         % solution initiale pour cas test �t�
uold(door) = dt;                      % prise en compte de la porte
uold(window) = wt;                    % prise en compte de la fen�tre

for l = 1:1000             % boucle en temps 
    %u = (eye(length(uold),length(uold))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
    u = (eye(length(uold),length(uold))+nu*deltat*A)*uold - nu*deltat*b; % explicite
    uold = u;
    
    pause(0.01)
    
    U = G;
    U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice
    % correspondant a la grille
    mesh(X,Y,U);                          % on dessine la solution sur la grille
    axis('ij');
end
%}