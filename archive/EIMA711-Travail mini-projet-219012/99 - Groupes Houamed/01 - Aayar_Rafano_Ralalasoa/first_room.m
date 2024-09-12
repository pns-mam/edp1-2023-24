%AAYAR - RAFANO - RALALASOA AVEC L'AIDE DE MME.DOLEAN
%%SCRIPT qui calcule la température ambiante dans une chambre
%%PARTIES 1 ET 2 INCLUSES 

n=30;         %nombre de points de la grille
wt = -10;     %température extérieure (celle près des fenêtres)
dt = 15;      % température de la porte
ht = 500;     % température du chauffage

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-0.9) & (X<-0.4) & (Y>-0.4) & (Y<0.9)) | ...  % geometrie de la chambre
   ((X>-0.4) & (X<0.2) & (Y>-1) & (Y<0.9)) | ... 
   ((X>0.2) & (X<0.9) & (Y>0.2) & (Y<0.9));
H = ((X>-0.9) & (X<-0.85) & (Y>0) & (Y<0.4))  %position du chauffage
%P = ((X>0.8) & (X<0.9) & (Y>0.4) & (Y<0.7))  %position de la porte 
%F = ((X>-0.3) & (X<0.1) & (Y>-1.1) & (Y<-0.9))   %position de la fenêtre
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  

%spy(G,"b");hold on;spy(H,"r");hold on; spy(P,"y");hold on;spy(F,"g") %pour vérifier la position 


A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
window = [];                        %indices correspondant a la fenetre
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              % CL Neumann sur le mur gauche
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i,j+1) == 0
        if (Y(i,j)>0.4) && (Y(i,j)<0.7) % On garde la CL Dirichlet pour la porte
          door = [door ind];
          else
        A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur droit
        end
      end
      if G(i+1,j) == 0 && i<n-1     
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas
      end
      if G(i-1,j) == 0 
        if (X(i,j)>-0.3 && X(i,j)<0.1)  % On garde la CL Dirichlet pour la fenêtre
          window = [window ind];
          else
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
          end
      end 
    end
  end
end

%%%%%%%%%%%%%%%%%%%PARTIE 1%%%%%%%%%%%%%%%%%%%%%%%

h = 2/(n-1);                         % pas du maillage
A = -A/h^2;                         % division par h^2
b = zeros(length(k),1);             %initialisation du second membre
heat = G(H);                        % indices correspondant au radiateur
b(window) = -1/h^2*wt;              % CL Dirichlet : fenetre
b(heat) = -ht;                      % prise en compte du chauffage
b(door) = -1/h^2*dt;                % CL Dirichlet : fenetre
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');     



%%%%%%%%%%%%%%%%%PARTIE 2%%%%%%%%%%%%%%%%%%
%{
h = 2/(n-1);                          % pas d'espace
A = -A/h^2;                           % division par h^2


alpha = 0.5; nu =1;                   % CFL 
deltat = alpha/(2*nu)*h^2;            
b = zeros(length(k),1);               % initialisation du second membre
heat = G(H);                          % indices correspondant au radiateur
b(heat) = -ht;                        % prise en compte du chauffage
b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*wt;                % CL Dirichet: fenetre
uold = zeros(length(k),1);            %solution initiale pour cas test ski
%uold = 35*ones(length(k),1);         % solution initiale pour cas test été
uold(door) = dt;                      % prise en compte de la porte
uold(window) = wt;                    % prise en compte de la fenêtre

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