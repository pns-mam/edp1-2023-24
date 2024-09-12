clear all;
close all;
%Chambre2SimuInstationnaire permet d'observer l'evolution de la temperature
%dans la chambre 2 avec un chauffage
%ot : temperature exterieure
%dt : temperature de la porte
%ht : temperature du chauffage
%n : nombre de points de la grille
function P2 = Chambre2SimuInstationnaire(ot,dt,ht,n)
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                            % grille bi-dimensionnelle 
G = ((Y>-0.75) & (Y<0.75) & (X>-0.9) & (X<0.5)) | ...  % geometrie de la chambre 2
   ((Y>0) & (Y<0.9) & (X>0.5) & (X<0.9)) | ...
   ((Y>0.75) & (Y<0.9) & (X>-0.9) & (X<-0.5)) | ... 
   ((Y>-0.9) & (Y<-0.75) & (X>0) & (X<0.5)) |...
   ((Y>0.1) & (Y<0.4) & (X>0.9) & (X<1));
%H = ((X>0.1) & (X<0.4) & (Y>-0.9) & (Y<-0.6));       % chauffagepres de la porte
H = ((X>-0.15) & (X<0.15) & (Y>-0.15) & (Y<0.15));      % chauffage au centre
%H = ((X>0.6) & (X<0.9) & (Y>0.1) & (Y<0.4));;      % chauffage pres de la fenetre
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0  
         A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur gauche
      end
      if G(i,j+1) == 0 && j<n-1         
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur droit
      end
      if G(i+1,j) == 0      
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
      end
      if G(i-1,j) == 0
        if (X(i,j)>0.1 && X(i,j)<0.4) % On garde la CL Dirichlet pour la porte 
          door = [door ind];  
        else
          A(ind,ind) = A(ind,ind)-1;   % CL Neumann sur le mur du haut
        end
      end 
    end
  end
end
h = 2/(n-1);
A = -A/h^2;                         % division par h^2
b = zeros(length(k),1);
window = G(G(:,end-1)>0,end-1);     % indices correspondant a la fenetre

heat = G(H);                        % indices correspondant au radiateur
b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(heat) = -ht;                      % et du chauffage
b(door) = -1/h^2*dt;

alpha = 0.5;
nu = 1;
deltat = alpha/(2*nu)*(h^2);
%u = 5*ones(length(k),1); %test pour une chambre froide initialement
u = 30*ones(length(k),1); %test pour une chambre chaude initialement

tic
for l = 1:5000
  uold = u;
  u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
  %u = (eye(size(A))+nu*deltat*A)*uold - nu*deltat*b; % explicite
  U = G;
  U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
  if mod(l,500) == 0
    mesh(X,Y,U);                        % correspondant a la grille 
    axis('ij');                         % on dessine la solution sur la grille
    pause(0.1);
  end
end                    
toc

c=0;
r=0;
for i=1:size(U,1)
  for j=1:size(U,2)
    if U(i,j)!=0
      r=r+U(i,j);
      c=c+1;
    end
  end
end
s=r/c;
disp(s);
disp(U);