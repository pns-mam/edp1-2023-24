# Mini projet EDP1
# Yoann Le Saint - Maxime Le

%clear all, close all, clc
% Simulation du probleme de Poisson dans le cas du probleme de Dirichlet homogene

n = 80;  % nombre de points de discretisation
wt = 20;                         % Temperature exterieur
dt = 20;                         % Temperature de la porte d'entree
ht = 150;                        % Temperature du chaffage

X=linspace(0,1,n);
[X,Y]=meshgrid(X,X);


G= ( ((X>0) & (X<0.9) & (Y>0) & (Y<0.6)) | ((X>0.9) & (X<1) & (Y>0.4) & (Y<0.6)) | ( (X>0.2) & (X<0.9) & (Y>0.6) & (Y<1) ) );

H = ((X>-0.80) & (X<0.91) & (Y>0.41) & (Y<0.59));  % Position du chauffage

k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A   

A = delsq(G);                       % Laplacien a l'interieur du domaine G
heat = G(H)>0;
door=[];
fenetre=[];

for i = 2:n-1                        
  for j = 2:n-1   
    ind = G(i,j);
    if G(i,j) ~= 0 
      if (G(i,j+1) == 0 || G(i+1,j) == 0)
        % mur de droite et mur du bas
        A(ind,ind)=A(ind,ind)-1;
      end
      if G(i,j-1) == 0
        % mur de gauche
        if (Y(i,j)>0.1 && Y(i,j)<0.5)
          fenetre = [fenetre ind];
        else A(ind,ind)=A(ind,ind)-1;
        end
      end
    if G(i-1,j) == 0                
        % mur du haut
        if (X(i,j)>0.24 && X(i,j)<0.31) 
            door = [door ind];
        else A(ind,ind)=A(ind,ind)-1;
        end
      end
    end
   end
end

h = 2/(n-1);
A = -A/h^2; 
b = zeros(length(k),1);


b(fenetre) = -(1/(h^2))*wt;              
b(heat) = -ht;                     
b(door) = -(1/(h^2))*dt;
u = A\b;                            
U = G;
U(G>0) = full(u(G(G>0)));



mesh(Y,X,U);
xlabel('Y');
ylabel('X');
axis('ij');     

moy=mean(mean(U))