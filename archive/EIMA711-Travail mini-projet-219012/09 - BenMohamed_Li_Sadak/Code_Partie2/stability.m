tic
%Simulation instationnaire
%Étude de la Stabilité 
%CHAMBRE NUMÉRO1
ht = 100; % temperature du chauffage
dt = 5;  % temperature de la porte d'entree
wt = -5; % temperature de la fenetre
n = 30;   % nombre de points de discretisation

x = linspace(-1,1,n);                      
[X,Y] = meshgrid(x,x);                  % grille bi-dimensionnelle 
G = ((X>0.3) & (X<0.8) & (Y>-0.7) & (Y<0.6)) | ...  % Partie de droite horizontal geometrie de la chambre
   ((X>0.35) & (X<0.75) & (Y<-0.7) & (Y>-0.8)) | ...   %fenetre haut gauche    
   ((X>-0.6) & (X<0.8) & (Y>-0.5) & (Y<0.6)) | ...  %Base de la chmbre
   ((X>0.3) & (X<0.8) & (Y<0.8) & (Y>0.6)) | ...   % base fenetre dans la base
   ((X>0.35) & (X<0.75) & (Y>0.8) & (Y<1)) | ...   % fenetre dans la base
   ((X<-0.6) & (X>-0.8) & (Y>-0.3) & (Y<0.4)) | ... %interieur de la base
   ((X>-0.9) & (X<-0.8) & (Y>-0.2) & (Y<0.35));    %porte
H = ((X>0.3) & (X<0.7) & (Y>-0.2) & (Y<0.35));
spy(G)

k = find(G);                          % les indices des G non nuls
G = zeros(size(G));                   % conversion logique - reel
G(k) = (1:length(k))';                % suite d'indices de la matrice A                           
heat = G(H); 

A = delsq(G);                         % matrice du -Laplacien sans 1/h^2 a l'interieur

door = [];                            % indices correspondant a la porte
window = [];                          % indices correspondant a la fenetre
for i = 2:n-1                        
  for j = 2:n-1
    ind = G(i,j);
    if ind ~= 0 
      if G(i-1,j) == 0                % la fenetre se situe sur le mur en haut
        if (X(i,j)>0.35 && X(i,j)<0.75) 
          window = [window G(i,j)];         % indices correspodant a la fenetre 
        else
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la porte le mur est isolant
        end
      end
      if G(i+1,j) == 0                % mur du bas isolant: Neumann
        if (X(i,j)>0.35 && X(i,j)<0.75)  
          window = [window G(i,j)];         % indices correspodant a la fenetre 
        else
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la porte le mur est isolant
        end
      end
      if G(i,j-1) == 0                % mur a gauche isolant: Neumann
         if (Y(i,j)>-0.2 && Y(i,j)<0.35) 
          door = [door ind];   % indices correspondant a la porte
        else
        A(ind,ind) = A(ind,ind)-1;
        end
      end
      if G(i,j+1) == 0                % la porte se situe sur le mur de droite 
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la fenetre le mur est isolant
      end
    end
  end
end

h=2/(n-1);
A = -A/h^2;   
b=zeros(length(k),1);
b(heat) = -ht;                        % prise en compte du chauffage
b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*wt;                % CL Dirichet: fenetre

%Étude de la stabilité schema Euler Explicite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stabilité avec 0<alpha<0.5
alpha=0.5;
%alpha=0.25; 
%alpha=0.1;  

%instabilité sinon
%alpha=0.6; 
%alpha=0.7;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stabilité pour tout alpha pour Euler Implicite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



nu=1;
deltat=alpha/(2*nu)*h^2;
%u=zeros(length(k),1);
u=5*ones(length(k),1);

for l=1:500
  uold=u;
  u=(eye(size(A))+nu*deltat*A)*uold - nu*deltat*b;
  %u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
  U=G;
  U(G>0)=full(u(G(G>0)));
  if mod(l,100)==0
    disp(l);
    mesh(X,Y,U);
    axis('ij');
    pause(0.1);
    end
end

toc
%disp(U)