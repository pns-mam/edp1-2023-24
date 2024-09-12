clear all

function chambre_h(ot,dt,ht,n)
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);    
                         % grille bi-dimensionnelle 
G = ((X>-0.9) & (X<-0.4) & (Y>-1) & (Y<0.8)) | ...  % geometrie de la chambre 
    ((X>-0.4) & (X<0.9) & (Y>-0.3) & (Y<0.3))| ...
    ((X>0.5) & (X<0.9) & (Y>0.3) & (Y<0.8))|...
    ((X>-0.8) & (X<-0.5) & (Y>0.8) & (Y<1))| ...
    ((X>0.6) & (X<0.8) & (Y>0.8) & (Y<1)); 
%H = ((X>-0.7) & (X<-0.6) & (Y>-0.2) & (Y<0.2))| ...       % 2 chauffages au centre
%    ((X>0.6) & (X<0.7) & (Y>-0.2) & (Y<0.2));      
H = ((X>-0.7) & (X<-0.6) & (Y>0.4) & (Y<0.8))| ...       % 2 chauffages en bas 
    ((X>0.6) & (X<0.7) & (Y>0.4) & (Y<0.8));   
  
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1 :length(k))';             % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0              % CL Neumann sur le mur gauche
        if (Y(i,j)>-0.2 && Y(i,j)<0.2) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        else
          A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i,j+1) == 0              % CL Neumann sur le mur droit
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0  && i<n-1    % On garde la CL Dirichlet pour la fenetre
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
endfunction


%realisation de nos cas tests stationnaire
%possibilite de faire varier la position du chauffage H en decommantant ci-dessus et en utilisant le cas 3
%ot temperature exterieur
%ht temperature chauffage
%dt temperature porte
%n maillage

%cas 1 : très chaud dehors et sans chauffage 
%chambre_h(20,20,0,30);

%cas 2 : très froid dehors et sans chauffage
chambre_h(-10,15,0,30);

%cas 3 : très froid dehors et avec chauffage
%chambre_h(-10,15,400,30);

%cas 4 : cas 3 avec un très grand n
%chambre_h(-10,15,400,400);









