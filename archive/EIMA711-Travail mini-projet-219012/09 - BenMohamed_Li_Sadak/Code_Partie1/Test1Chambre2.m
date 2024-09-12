% Simulation du probleme de Poisson
%Température ambiante en été avec les portes et fenêtres à 20°
%CHAMBRE NUMÉRO 2
ht = 0; % pas de temperature du chauffage
dt = 20;  % temperature de la porte d'entree à 20°
wt = 20; % temperature de la fenetre à 20°
n = 30;   % nombre de points de discretisation

x = linspace(-1,1,n);                      
[X,Y] = meshgrid(x,x);       
G = ((X>-0.8) & (X<0.8) & (Y>-0.5) & (Y<0.5)) | ...  % grille bi-dimensionnelle 
   ((X>-0.5) & (X<0.5) & (Y>-0.8) & (Y<0.8))| ...
  ((X>-0.9) & (X<-0.8) & (Y>-0.2) & (Y<0.2)) | ...
  ((X<0.9) & (X>0.8) & (Y>-0.2) & (Y<0.2)) | ...
  ((X>-0.2) & (X<0.2) & (Y>0.8) & (Y<1));
H = ((X>-0.2) & (X<0.2) & (Y>-0.6) & (Y<-0.4));  %position du chauffage
spy(G)

 
k = find(G);                          % les indices des G non nuls
G = zeros(size(G));                   % conversion logique - reel
G(k) = (1:length(k))';                % suite d'indices de la matrice A
spy(G)                                % on dessine la grille
heat = G(H); 
A = delsq(G);                         % matrice du -Laplacien sans 1/h^2 a l'interieur
door = [];                            % indices correspondant a la porte
window = [];                          % indices correspondant a la fenetre

for i = 2:n-1                        
  for j = 2:n-1
    ind = G(i,j);
    if ind ~= 0 
      if G(i-1,j) == 0                % la porte se situe sur le mur du haut
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la porte le mur est isolant
      end
      if G(i+1,j) == 0                % mur du bas isolant: Neumann
        if (X(i,j)>-0.2 && X(i,j)<0.2) 
          door = [door ind];   % indices correspodant a la fenetre 
        else
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la fenetre le mur est isolant
        end
      end
      if G(i,j-1) == 0                % mur a gauche isolant: Neumann
        if (Y(i,j)>-0.2 && Y(i,j)<0.2) 
          window = [window ind];   % indices correspodant a la fenetre 
        else
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la fenetre le mur est isolant
        end
      end
      if G(i,j+1) == 0                % la fenetre se situe sur le mur a droite
        if (Y(i,j)>-0.2 && Y(i,j)<0.2) 
          window = [window ind];   % indices correspodant a la fenetre 
        else
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la fenetre le mur est isolant
        end
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