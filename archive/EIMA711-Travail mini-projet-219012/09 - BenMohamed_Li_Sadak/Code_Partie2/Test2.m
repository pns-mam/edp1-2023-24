tic
%Simulation instationnaire numéro 2
%Temperature extérieure très élevée et on met le radiateur en mode ”clim”.
% Combien de temps il faut pour atteindre la temperature id́eale (17°-22°)?
%CHAMBRE NUMÉRO1
ht = -300; % temperature du chauffage en mode clim pour refroidir
dt = 29;  % temperature de la porte d'entree très élevée
wt = 35; % temperature de la fenetre très élevée
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

alpha=0.4;
nu=1;
deltat=alpha/(2*nu)*h^2;
%u=zeros(length(k),1);
u=30*ones(length(k),1);               %temperature initiale très élevée


for l=1:1000
  uold=u;
  u=(eye(size(A))+nu*deltat*A)*uold - nu*deltat*b;
  %u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b); % implicite
  U=G;
  U(G>0)=full(u(G(G>0)));
  %mesh(X,Y,U);
  moyenne=sum(U(:))/(length(find(U)));
  if mod(l,100)==0
    mesh(X,Y,U);
    axis('ij');
    %disp(l);  %nb d'iterations
    disp("temperature");
    disp(moyenne);
    pause(1);
    %pause(10); % pour voir plus doucement les differentes états qui sont dans le rapport
  end
    if 17<moyenne & moyenne<22
      mesh(X,Y,U);
      axis('ij');
      disp("temperature idéal atteinte")
      disp(moyenne);
      break;
    end
end
toc
%disp(U)