function [A,G,b,door,window] = Param_partie2(num_chambre, ot,dt,ht,n)
  x = linspace(-1,1,n);                            
  [X,Y] = meshgrid(x,x);  
if( num_chambre == 1)  
  G1 = ((X>-1) & (X<1) & (Y>-0.33) & (Y<0.8)) | ((X>-1) & (X<-0.5) & (Y>-1) & (Y<-0.33)) | ((X>-0.6) & (X<0) & (Y>0.8) & (Y<1));
  H = ((X>-1) & (X<-0.7) & (Y>-0.3) & (Y<0.3));     
  k = find(G1);                        
  G = zeros(size(G1));                 
 
 
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
        if (X(i,j)>-0.9 && X(i,j)<-0.7) 
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
end
if(num_chambre == 2)
G1 = ((X>-0.7) & (X<0.7) & (Y>0) & (Y<0.9)) | ((X>-0.5) & (X<0.5) & (Y>-0.8) & (Y<0)) ;

H = ((X>-0.2) & (X<0.2) & (Y>-0.5) & (Y<-0.2));
k = find(G1);                        % les indices des G non nuls
G = zeros(size(G1));                 % conversion logique - reel 
G(k) = (1:length(k))';              % suite d'indices de la matrice A  

A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = [];                          % indices correspondant a la porte
window = [];
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i+1,j) == 0              % mur du haut isolant donc Neumann
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i-1,j) == 0              % mur du bas
        if (X(i,j)>-0.25 && X(i,j) < 0.25)
          window = [window ind]; %fenetre
        else
          A(ind,ind) = A(ind,ind)-1; %sinon isolant
        end
        
      end
      if G(i,j-1) == 0  %mur a gauche
        if (Y(i,j) >0.3 && Y(i,j) < 0.8)
          door = [door ind];
        else
          A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
        end
        
      end
      if G(i,j+1) == 0 %mur a droite
          A(ind,ind) = A(ind,ind)-1;
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


end
