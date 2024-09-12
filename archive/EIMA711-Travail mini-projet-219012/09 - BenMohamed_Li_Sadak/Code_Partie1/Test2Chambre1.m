% Simulation du probleme de Poisson
%empérature ambiante en hiver, sans chauffage, -10°à l’extérieur, et lesportes à 15°
%CHAMBRE NUMÉRO1
ht = 0; % pas de temperature du chauffage
dt = 15;  % temperature de la porte d'entree à 20°
wt = -10; % temperature de la fenetre à 20°
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
H = ((X>0.3) & (X<0.7) & (Y>-0.2) & (Y<0.35)); %position du chauffage
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


h = 2/(n-1);                          % pas du maillage
A = -A/h^2;                           % division par h^2
b = zeros(length(k),1);               % initialisation du second membre
b(heat) = -ht;                        % prise en compte du chauffage
b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*wt;                % CL Dirichet: fenetre
u = A\b;                              % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice 
mesh(X,Y,U);                          % on dessine la solution sur la grille
axis('ij');

disp("Temperature moyenne de la chambre")
moyenne=sum(U(:))/(length(find(U)));
disp(moyenne)