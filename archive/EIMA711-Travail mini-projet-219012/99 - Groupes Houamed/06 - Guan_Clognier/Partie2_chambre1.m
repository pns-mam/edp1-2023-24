clear all, close all, clc
ot=-20;
dt=20;
ht=0;
n=40;
% ROOMTEMPERATURE calcule la température ambiante dans une chambre
% U = RoomTemperature(ot, dt, ht, n); prend la température extérieure ot,
% la température de la porte dt et la température du chauffage ht et la
% nombre de points de grille n, puis calcule la température ambiante de la
% chambre

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X>-0.8) & (X<0.2) & (Y>-0.6) & (Y<0.6)) | ...  % geometrie de la chambre
   ((Y>(-0.2)*X-0.76)&(Y<0.2*X+0.76)&(X>-0.8) & (X<0.2))|...
   ((X>0.2) & (X<0.6) & (Y>-0.5) & (Y<0.6)) | ... 
   ((X>-1) & (X<-0.8) & (Y>-0.5) & (Y<-0.3)) | ... 
   ((X>-1) & (X<-0.8) & (Y>0.4) & (Y<0.6));
%H = ((X>-0.6) & (X<0) & (Y>0.5) & (Y<0.75));       % position du chauffage
H = ((X>-0.45) & (X<-0.05) & (Y>-0.2) & (Y<0.2));
%H = ((X>0.4) & (X<0.6) & (Y>-0.3) & (Y<0.3));      % chauffage près le mur
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k));              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0  && j>4           % CL Neumann sur le mur gauche
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i,j+1) == 0              % CL Neumann sur le mur droit
        if (Y(i,j)>-0.2 && Y(i,j)<0.3) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        else
          A(ind,ind) = A(ind,ind)-1;
         end
      end
      if G(i+1,j) == 0     % On garde la CL Dirichlet pour la fenetre
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
      end
      if G(i-1,j) == 0
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut 
      end 
    end
  end
end
h = 2/(n-1);
A = -A/h^2;                         % division par h^2
b = zeros(length(k),1);
window = G(G(:,2)>0,2);     % indices correspondant a la fenetre
heat = G(H);                        % indices correspondant au radiateur
b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(heat) = -ht;                      % et du chauffage
b(door) = -1/h^2*dt;


alpha = 0.5; nu =1;
deltat = alpha/(2*nu)*h^2;
u = A\b;
b(heat)=-500;
tic
for l = 1:15000               % boucle en temps 
    uold = u;
    u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
    %u = (eye(size(A))+nu*deltat*A)*uold - nu*deltat*b; % explicite
    U = G;
    U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice
    % correspondant a la grille
    if(mod(l,100)==0)
    time=deltat*l;
    mesh(X,Y,U);                          % on dessine la solution sur la grille
    title('time=',num2str(time));
    axis('ij');
    pause(0.1);
    end
end
mytimer1=toc;
disp(mytimer1);
%mesh(X,Y,T);                        % correspondant a la grille 
%axis('ij');                         % on dessine la solution sur la grille
