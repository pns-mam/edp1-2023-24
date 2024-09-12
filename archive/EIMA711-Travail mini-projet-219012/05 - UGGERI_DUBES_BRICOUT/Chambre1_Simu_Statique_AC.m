clear all;
close all;
%Chambre1SimuStatiqueAC permet de calculer la temperature ambiante dans la chambre 1 avec un chauffage
%ot : temperature exterieure
%dt : temperature de la porte
%ht : temperature du chauffage
%n : nombre de points de la grille
function P1 = Chambre1SimuStatiqueAC(ot,dt,ht,n)
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                           % grille bi-dimensionnelle 
G = (((X>-0.9) & (X<0.9) & (Y>-0.9) & (Y<0))|... % geometrie de la chambre 1
     ((X>-0.2) & (X<0.9) & (Y>-0.8) & (Y<0.9))|...
     ((X>0.9) & (X<1) & (Y>-0.9) & (Y<-0.5))) ;
H = ((X>0) & (X<0.5) & (Y>-0.5) & (Y<0));       % chauffage au centre de la chambre
%H = ((X>0.5) & (X<0.9) & (Y>-0.9) & (Y<-0.5));      % chauffage pres de la fenetre
%H = ((X>0.1) & (X<0.4) & (Y>-0.2) & (Y<0.2));      % chauffage pres de la porte
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
        if (Y(i,j)>0.1 && Y(i,j)<0.4)  
          door = [door ind];      % On garde la CL Dirichlet pour la porte
        else
          A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur gauche
        end
      end
      if G(i,j+1) == 0 && j<n-1         
        A(ind,ind) = A(ind,ind)-1;    % CL Neumann sur le mur droit
      end
      if G(i+1,j) == 0      
        A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
      end
      if G(i-1,j) == 0
        A(ind,ind) = A(ind,ind)-1;       % CL Neumann sur le mur du haut
         
      end 
    end
  end
end
h = 2/(n-1);
A = -A/h^2;                         % division par h^2
b = zeros(length(k),1);
window = G(G(:,end-1)>0,end-1);     % indices correspondant a la fenetre

heat = G(H);                        % indices correspondant au radiateur faire attention qu'il se trouve bien dans la piece
b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet pour la fenetre
b(heat) = -ht;                      % prise en compte du chauffage
b(door) = -1/h^2*dt;
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij');  


%boucle qui permet de calculer la temperature ambiante d'une piece
%en faisant la moyenne des temperatures des differents endroits de la piece
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