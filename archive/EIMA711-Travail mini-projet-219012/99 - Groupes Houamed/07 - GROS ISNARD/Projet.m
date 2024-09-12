clear, clc

ht = 500; % temperature du chauffage
dt = 15;  % temperature de la porte d'entree
wt = -10; % temperature de la fenetre
ot = -10;  % temperature dehors
n = 50;  % nombre de points de discretisation

cfl = 0.25;  %cfl max = 0.25 en explicite
nt = 160;

%[U,u]=RoomTemperatureStationnaire1(ot,dt,ht,n)
%[U,u]=RoomTemperatureStationnaire2(ot,dt,ht,n)
%[Ue]=RoomTemperatureExplicite1(ot,dt,ht,n,nt,cfl)  % en mode chauffage
%[Ue]=RoomTemperatureExplicite2(ot,dt,ht,n,nt,cfl)  % en mode clim
%[Ui]=RoomTemperatureImplicite1(ot,dt,ht,n,nt,cfl)  % en mode chauffage
%[Ui]=RoomTemperatureImplicite2(ot,dt,ht,n,nt,cfl)  % en mode clim

function [U,u] = RoomTemperatureStationnaire1(ot,dt,ht,n)

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-1) & (X<0) & (Y>-1) & (Y<0.95)) | ...  % geometrie de la chambre
   ((X>0) & (X<1) & (Y>0) & (Y<0.95)) | ... 
   ((X>0.3) & (X<0.7) & (Y>-0.2) & (Y<0)) | ...
   ((X>-0.7) & (X<-0.2) & (Y>0.95) & (Y<1)) | ...
   ((X>0.2) & (X<0.7) & (Y>0.95) & (Y<1));
%H=[];  %pas de chauffage
H = ((X>-1) & (X<-0.8) & (Y>-0.1) & (Y<0.3));       % position du chauffage
%H =  ((X>-0.2) & (X<0.2) & (Y>0.75) & (Y<0.95));    %Meilleure position
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
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
      if G(i+1,j) == 0 && i<n-1     % On garde la CL Dirichlet pour la fenetre
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
      end
      if G(i-1,j) == 0
        if (X(i,j)>0.3 && X(i,j)<0.7 && Y(i,j)>-0.2 && Y(i,j)<-0.1 ) 
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
window = G(end-1,G(end-1,:)>0);     % indices correspondant a la fenetre
heat = G(H);                        % indices correspondant au radiateur
b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(heat) = -ht;                      % et du chauffage
b(door) = -1/h^2*dt;
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');                         % on dessine la solution sur la grille
end


function [U,u] = RoomTemperatureStationnaire2(ot,dt,ht,n)
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                            
G = ((Y<0.8) & (Y<1) & (X>-1) & (X<1) & (Y>-1)) | ((Y>0.8) & (Y<1) & (X >-0.2) & (X<0.2)); % geometry of the room
H = ((Y>-0.8) & (Y<-0.4) & (X<-0.8) & (X>-1)); % heat

k = find(G);                       
G = zeros(size(G));                 
G(k) = (1:length(k))';
A = delsq(G);
door = [];                          
window = [];                        
for i = 2:n-1                       
  for j = 2:n-1                     
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              
        if ((Y(i,j)>0.2) && (Y(i,j)<0.8))
          window = [window ind];
        else
            A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i,j+1) == 0
        A(ind,ind)=A(ind,ind)-1;
      end
      if G(i+1,j) == 0   
        if (X(i,j)>-0.2 && X(i,j)<0.2)
            door = [door ind];
        else
          A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i-1,j) == 0
        A(ind,ind)=A(ind,ind)-1;
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
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice

spy(G)
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');  
xlabel('X')
ylabel('Y')
zlabel('Temperature')
end

function [U] = RoomTemperatureExplicite1(ot,dt,ht,n,nt,cfl)

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-1) & (X<0) & (Y>-1) & (Y<0.95)) | ...  % geometrie de la chambre
   ((X>0) & (X<1) & (Y>0) & (Y<0.95)) | ... 
   ((X>0.3) & (X<0.7) & (Y>-0.2) & (Y<0)) | ...
   ((X>-0.7) & (X<-0.2) & (Y>0.95) & (Y<1)) | ...
   ((X>0.2) & (X<0.7) & (Y>0.95) & (Y<1));
%H=[];  %pas de chauffage
H = ((X>-1) & (X<-0.8) & (Y>-0.1) & (Y<0.3));       % position du chauffage
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
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
      if G(i+1,j) == 0 && i<n-1     % On garde la CL Dirichlet pour la fenetre
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
      end
      if G(i-1,j) == 0
        if (X(i,j)>0.3 && X(i,j)<0.7 && Y(i,j)>-0.2 && Y(i,j)<-0.1 ) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        else
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
        end
      end 
    end
  end
end

h = 2/(n-1);
t = h^2*cfl;

b = zeros(length(k),1);
window = G(end-1,G(end-1,:)>0);     % indices correspondant a la fenetre
heat = G(H);                        % indices correspondant au radiateur
b(window) = 1/h^2*ot*t;              % prise en compte des CL Dirichlet
b(heat) = ht*t;                      % et du chauffage
b(door) = 1/h^2*dt*t;

[U,u]=RoomTemperatureStationnaire1(ot,dt,0,n);

B= -A*(t/h^2) + eye(length(k));

for i=1:nt
  u = B*u + b;                            % solution du systeme sous forme de vecteur
  U = G;
  U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
  if mod(i,20)==0
    mesh(X,Y,U);                        % correspondant a la grille 
    axis('ij');                         % on dessine la solution sur la grille
    pause(1);
  end
end
end

function [U] = RoomTemperatureExplicite2(ot,dt,ht,n,nt,cfl)
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                            
G = ((Y<0.8) & (Y<1) & (X>-1) & (X<1) & (Y>-1)) | ((Y>0.8) & (Y<1) & (X >-0.2) & (X<0.2)); % geometry of the room
H = ((Y>-0.8) & (Y<-0.4) & (X<-0.8) & (X>-1)); % heat

k = find(G);                       
G = zeros(size(G));                 
G(k) = (1:length(k))';
A = delsq(G);
door = [];                          
window = [];                        
for i = 2:n-1                       
  for j = 2:n-1                     
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              
        if ((Y(i,j)>0.2) && (Y(i,j)<0.8))
          window = [window ind];
        else
            A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i,j+1) == 0
        A(ind,ind)=A(ind,ind)-1;
      end
      if G(i+1,j) == 0   
        if (X(i,j)>-0.2 && X(i,j)<0.2)
            door = [door ind];
        else
          A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i-1,j) == 0
        A(ind,ind)=A(ind,ind)-1;
      end 
    end
  end
end
h = 2/(n-1);
t = h^2*cfl;

b = zeros(length(k),1);
heat = G(H);                        % indices correspondant au radiateur
b(window) = 1/h^2*ot*t;              % prise en compte des CL Dirichlet
b(heat) = ht*t;                      % et du chauffage
b(door) = 1/h^2*dt*t;

init=ones(length(k),1)*30 ;     %pour la température homogène en été

u=init;

B= -A*(t/h^2) + eye(length(k));
for i=1:nt
  u = B*u + b;                            % solution du systeme sous forme de vecteur
  if mod(i,20)==0
    U = G;
    U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
    mesh(X,Y,U);                        % correspondant a la grille 
    axis('ij');                         % on dessine la solution sur la grille
    pause(1);
  end
end

end

function [U] = RoomTemperatureImplicite1(ot,dt,ht,n,nt,cfl)
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-1) & (X<0) & (Y>-1) & (Y<0.95)) | ...  % geometrie de la chambre
   ((X>0) & (X<1) & (Y>0) & (Y<0.95)) | ... 
   ((X>0.3) & (X<0.7) & (Y>-0.2) & (Y<0)) | ...
   ((X>-0.7) & (X<-0.2) & (Y>0.95) & (Y<1)) | ...
   ((X>0.2) & (X<0.7) & (Y>0.95) & (Y<1));
%H=[];  %pas de chauffage
H = ((X>-1) & (X<-0.8) & (Y>-0.1) & (Y<0.3));       % position du chauffage
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
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
      if G(i+1,j) == 0 && i<n-1     % On garde la CL Dirichlet pour la fenetre
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
      end
      if G(i-1,j) == 0
        if (X(i,j)>0.3 && X(i,j)<0.7 && Y(i,j)>-0.2 && Y(i,j)<-0.1 ) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        else
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
        end
      end 
    end
  end
end

h = 2/(n-1);
t = h^2*cfl;

b = zeros(length(k),1);
window = G(end-1,G(end-1,:)>0);     % indices correspondant a la fenetre
heat = G(H);                        % indices correspondant au radiateur
b(window) = 1/h^2*ot*t;              % prise en compte des CL Dirichlet
b(heat) = ht*t;                      % et du chauffage
b(door) = 1/h^2*dt*t;

[U,u]=RoomTemperatureStationnaire1(ot,dt,0,n);

B= A*(t/h^2) + eye(length(k));
L=chol(B);                           %car matrice symétrique défini positive
for i=1:nt
  w = L'\(u+b);
  u = L\w;
  U = G;
  U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
  if mod(i,20)==0
    mesh(X,Y,U);                        % correspondant a la grille 
    axis('ij');                         % on dessine la solution sur la grille
    pause(1);
  end
end
end

function [U] = RoomTemperatureImplicite2(ot,dt,ht,n,nt,cfl)
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                            
G = ((Y<0.8) & (Y<1) & (X>-1) & (X<1) & (Y>-1)) | ((Y>0.8) & (Y<1) & (X >-0.2) & (X<0.2)); % geometry of the room
H = ((Y>-0.8) & (Y<-0.4) & (X<-0.8) & (X>-1)); % heat

k = find(G);                       
G = zeros(size(G));                 
G(k) = (1:length(k))';
A = delsq(G);
door = [];                          
window = [];                        
for i = 2:n-1                       
  for j = 2:n-1                     
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              
        if ((Y(i,j)>0.2) && (Y(i,j)<0.8))
          window = [window ind];
        else
            A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i,j+1) == 0
        A(ind,ind)=A(ind,ind)-1;
      end
      if G(i+1,j) == 0   
        if (X(i,j)>-0.2 && X(i,j)<0.2)
            door = [door ind];
        else
          A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i-1,j) == 0
        A(ind,ind)=A(ind,ind)-1;
      end 
    end
  end
end

h = 2/(n-1);
t = h^2*cfl;

b = zeros(length(k),1);
heat = G(H);                        % indices correspondant au radiateur
b(window) = 1/h^2*ot*t;              % prise en compte des CL Dirichlet
b(heat) = ht*t;                      % et du chauffage
b(door) = 1/h^2*dt*t;

init=ones(length(k),1)*30 ;     %pour la température homogène en été

u=init;

B= A*(t/h^2) + eye(length(k));
L=chol(B);
for i=1:nt
  w = L'\(u+b);
  u = L\w;
  if mod(i,20)==0
    U = G;
    U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
    mesh(X,Y,U);                        % correspondant a la grille 
    axis('ij');                         % on dessine la solution sur la grille
    pause(1);
  end
  
end

end

