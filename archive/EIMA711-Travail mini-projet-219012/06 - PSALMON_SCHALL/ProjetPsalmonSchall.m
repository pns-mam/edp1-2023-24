clear
clearvars

function D = delsq(G)
[m,n] = size(G);

% Indices of interior points
p = find(G);

% Connect interior points to themselves with 4's.
i = G(p);
j = G(p);
s = 4*ones(size(p));

% for k = north, east, south, west
for k = [-1 m 1 -m]
   % Possible neighbors in k-th direction
   Q = G(p+k);
   % Index of points with interior neighbors
   q = find(Q);
   % Connect interior points to neighbors with -1's.
   i = [i; G(p(q))];
   j = [j; Q(q)];
   s = [s; -ones(length(q),1)];
end
D = sparse(i,j,s);
endfunction

%Calcul de la solution stationnaire pour la chambre 1
function U = RoomTemperature1(ot,dt,ht,n)
  
  
  x = linspace(-1,1,n);                            
  [X,Y] = meshgrid(x,x);                                % grille bi-dimensionnelle 
  G = ((X>-0.9) & (X<0.9) & (Y>-0.3) & (Y<1) | ...      % geometrie de la chambre
      (X>0.1) & (X<0.9) & (Y>-1) & (Y<-0.3));
  H = ((X>-0.85) & (X<0.05) & (Y>-0.3) & (Y<-0.1));     % position du chauffage
  
  
  
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
        if G(i,j-1) == 0              % CL Neumann sur le mur gauche
          A(ind,ind) = A(ind,ind)-1;
        end
        if G(i,j+1) == 0              % CL Neumann sur le mur droit
          A(ind,ind) = A(ind,ind)-1;
        end
        if G(i+1,j) == 0     
          if (X(i,j)>-0.2 && X(i,j)<0.7) 
            window = [window ind];        % On garde la CL Dirichlet pour la fenetre 
          else
            A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
          end
        end
        if G(i-1,j) == 0
          if (X(i,j)>0.3 && X(i,j)<0.7) 
            door = [door ind];        % On garde la CL Dirichlet pour la porte 
          else
            A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
          end
        end 
      end
    end
  end
  h = 2/(n-1);
  A = -A/h^2;                      % division par h^2
  b = zeros(length(k),1);
  heat = G(H);                        % indices correspondant au radiateur
  b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
  b(heat) = -ht;                      % et du chauffage
  b(door) = -1/h^2*dt;
  u = A\b;                            % solution du systeme sous forme de vecteur
  U = G;
  U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
  subplot(2,2,1);                     % on met le graph en haut a gauche
  mesh(X,Y,U);                        % correspondant a la grille 
  axis('ij');                         % on dessine la solution sur la grille

  endfunction

%Calcul de la solution stationnaire pour la chambre 2
function U = RoomTemperature2(ot,dt,ht,n)

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-0.8) & (X<-0.2) & (Y>-0.8) & (Y<0.5) | ...  % geometrie de la chambre
    ((X>-0.2) & (X<0.8) & (Y>-0.7) & (Y<0.9)));    

H= ((X>0.65) & (X<0.8) & (Y>-0.1) & (Y<0.3) |...
((X>-0.7) & (X<-0.4) & (Y>-0.3) & (Y<0)));         % position des chauffages


k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G

door = [];                          % indices correspondant a la porte
window=[];                          % indices correspondant a la fenetre

for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              
        if (Y(i,j)>-0.3 && Y(i,j)<0)
          window=[window ind];      % On garde la CL Dirichlet pour la fenetre
        else
          A(ind,ind) = A(ind,ind)-1;% CL Neumann sur le mur gauche
        end 
    end    
      if G(i,j+1) == 0              % CL Neumann sur le mur droit
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0      
            A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du bas
      end
      if G(i-1,j) == 0
        if (X(i,j)>0.2 && X(i,j)<0.5) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        
        else
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
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
b(heat) = -ht;                      % et des chauffages
b(door) = -1/h^2*dt;
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
subplot(2,2,2);                     % on met le graph en haut a droite
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');                         % on dessine la solution sur la grille

endfunction

%Calcul de la solution instationnaire pour la chambre 1
function U = RoomTemperatureInstationaire1(ot,dt,ht,n)
  x = linspace(-1,1,n);                            
  [X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
  G = ((X>-0.9) & (X<0.9) & (Y>-0.3) & (Y<1) | ...   % geometrie de la chambre
      (X>0.1) & (X<0.9) & (Y>-1) & (Y<-0.3));
      
  H = ((X>-0.85) & (X<0.05) & (Y>-0.3) & (Y<-0.1));  % position du chauffage
  
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
        if G(i,j-1) == 0              % CL Neumann sur le mur gauche
          A(ind,ind) = A(ind,ind)-1;
        end
        if G(i,j+1) == 0              % CL Neumann sur le mur droit
          A(ind,ind) = A(ind,ind)-1;
        end
        if G(i+1,j) == 0     
          if (X(i,j)>-0.2 && X(i,j)<0.7) 
            window = [window ind];        % On garde la CL Dirichlet pour la fenetre 
          else
            A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
          end
        end
        if G(i-1,j) == 0
          if (X(i,j)>0.3 && X(i,j)<0.7) 
            door = [door ind];        % On garde la CL Dirichlet pour la porte 
          else
            A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
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
  alpha = 10; nu =1;
  deltat = alpha/(2*nu)*h^2;     % calcul du delta temps 
  u = ones(length(k),1)*-10;     % solution initiale
  for l = 1:150                  % boucle en temps 
      uold = u;
      u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
      %u = (eye(size(A))+nu*deltat*A)*uold - nu*deltat*b;  % explicite
      U = G;
      U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice
                                            % correspondant a la grille
      subplot(2,2,3);                       % on met le graph en bas a gauche
      mesh(X,Y,U);                          % on dessine la solution sur la grille
      axis('ij');
      pause(0.001)
  end

endfunction

%Calcul de la solution instationnaire pour la chambre 2
function U = RoomTemperatureInstationaire2(ot,dt,ht,n)
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-0.8) & (X<-0.2) & (Y>-0.8) & (Y<0.5) | ...  % geometrie de la chambre
    ((X>-0.2) & (X<0.8) & (Y>-0.7) & (Y<0.9)));    

H= ((X>0.65) & (X<0.8) & (Y>-0.1) & (Y<0.3) |...
((X>-0.7) & (X<-0.4) & (Y>-0.3) & (Y<0)));         % position des chauffages


k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G

door = [];                          % indices correspondant a la porte
window=[];                          % indices correspondant a la fenetre
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              
        if (Y(i,j)>-0.3 && Y(i,j)<0)
          window=[window ind];      % On garde la CL Dirichlet pour la fenetre
        else
          A(ind,ind) = A(ind,ind)-1;% CL Neumann sur le mur gauche
        end 
    end    
      if G(i,j+1) == 0              % CL Neumann sur le mur droit
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0      
            A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du bas
      end
      if G(i-1,j) == 0
        if (X(i,j)>0.2 && X(i,j)<0.5) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        
        else
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
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
  b(heat)= -ht;                       % et des chauffages
  b(door) = -1/h^2*dt;
  alpha = 10; nu =1;
  deltat = alpha/(2*nu)*h^2;     % calcul du delta temps 
  u = ones(length(k),1)*-10;     % solution initiale
  for l = 1:150                  % boucle en temps 
      uold = u;
      u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
      %u = (eye(size(A))+nu*deltat*A)*uold - nu*deltat*b; % explicite
      U = G;
      U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice
                                            % correspondant a la grille
      subplot(2,2,4);                       % on met le graph en bas a droite 
      mesh(X,Y,U);                          % on dessine la solution sur la grille
      axis('ij');
      pause(0.001)
  end

endfunction

RoomTemperature1(-10,15,180,30);
RoomTemperature2(-10,15,180,30);
RoomTemperatureInstationaire1(-10,15,180,30);
RoomTemperatureInstationaire2(-10,15,180,30);



