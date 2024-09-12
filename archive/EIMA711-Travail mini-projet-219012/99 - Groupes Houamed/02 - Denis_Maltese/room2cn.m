function U = room2cn(ot,dt,ht,hdt,n,nu,alpha,eps,tinit)
  %{
Calcul de la solution it�rative du probl�me de la partie 2 
pour notre chambre 2
en utilisant la m�thode de Crank-Nicolson avec:
ot:la temp�rature ext�rieure(en �C)
dt:la temp�rature � la porte(en �C)
ht et hdt: les "puissances " des deux chauffages(sans unit�)
n: le nombre pas de discr�tisation en espace pour la grille 2d
dans laquelle est dessin�e notre chambre 
nu: vitesse de diffusion(dx/dt car dx=dy)
alpha: nombre de cfl
eps: tol�rance de convergence, valeur de la norme de u-uold
en dessous de laquelle on consid�re que le mod�le a converg�
tinit:temp�rature uniforme dans la pi�ce avant le d�but de
la simulation
%}
h = 2/(n-1);   %longueur du pas d'espace
x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                          % grille bi-dimensionnelle 
G = ((X>-1) & (X<-0.2) & (Y>-1) & (Y<1)) | ...  % geometrie de la chambre
   ((X>0.3) & (X<1) & (Y>-1) & (Y<1)) | ... 
   ((X>-0.3) & (X<0.3) & (Y>0.4) & (Y<1));
H = ((X>-0.5) & (X<-0.2) & (Y>-0.8) & (Y<0.2));       % position du chauffage 1
H2 = ((X>0.2) & (X<0.5) & (Y>0.4) & (Y<0.7));         % position du chauffage 2
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G 
window1=[];              %indices correspondant aux fen�tres
window2=[];
window3=[]; 
door = [];               % indices correspondant a la porte
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0 % CL Neumann sur le mur gauche
        if (Y(i,j)>-0.7&& Y(i,j)<0 && X(i,j)<0)
          window1 = [window1 ind];   % On garde la CL Dirichlet pour la fen�tre
        else
          A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i,j+1) == 0              % CL Neumann sur le mur droit
        if (Y(i,j)>-0.9 && Y(i,j)<-0.3 && X(i,j)>0)
          window2 = [window2 ind];  % On garde la CL Dirichlet pour la fen�tre
        else
          A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i-1,j) == 0              % CL Neumann sur le mur du haut
        A(ind,ind) = A(ind,ind)-1;   
      end
      if G(i+1,j) == 0              % CL Neumann sur le mur du bas
        if (X(i,j)>-0.15 && X(i,j)<0.15) 
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        elseif (X(i,j)>0.35&& X(i,j)<0.95)
          window3 = [window3 ind];  % On garde la CL Dirichlet pour la fen�tre
        else
          A(ind,ind) = A(ind,ind)-1; 
        end
      end 
    end
  end
end
A = -A/h^2;                         % division par h^2
b = zeros(length(k),1); 
heat = G(H);                        % indices correspondant aux radiateurs
heat2=G(H2);                        
b(window1) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(window2) = -1/h^2*ot;   
b(window3) = -1/h^2*ot;   
b(door) = -1/h^2*dt;
b(heat) = -ht;                      % et des chauffages
b(heat2)=-hdt;
deltat = alpha/(2*nu)*h^2;           %pas en temps

u = tinit*ones(length(k),1);     % condition initiale
for l = 1:10000           % boucle en temps 
    uold = u;
    u = (eye(size(A))-((nu*deltat)/2)*A)\((eye(size(A))+((nu*deltat)/2)*A)*uold - nu*deltat*b); % explicite
    U = G;
    U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice
    mesh(X,Y,U);                          % on dessine la solution sur la grille
    axis('ij');
    if norm(u-uold)<eps                   %condition de convergence
      disp("convergence au pas de temps:")
      disp(l)
      return
    end
    pause(0.01)
end