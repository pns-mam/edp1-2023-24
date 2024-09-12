function U =Fonction_RoomTemperature2(ot,dt,ht,ht1,n)

% La fonction calcule la température ambiante dans la chambre 2 pour la partie 1 du projet

% Température fenêtre : ot (Température extérieure)
% Température porte : dt 
% Température du premier chauffage : ht
% Température du deuxième chauffage : ht1
% Nombre de points de la grille : n


% On dessine la chambre dans un premier temps avec les conditions de Neumann et de Dirichlet 
% puis on place le chauffage, la porte et la fenêtre
% on résout l'équation de poisson 
% Et on affiche la température ambiante de la chambre 2



x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             

G  = ((X>-1) & (X<0.3) & (Y>-0.4) & (Y<0.6)) | ((X>0.3) & (X<1) & (Y>-1) & (Y<1)); % Chambre 2

H= ((X>0.8) & (X<1) & (Y>-0.2) & (Y<0.2))           % positon du premier chauffage (ht)
H1 = ((X>-0.7) & (X<-0.5) & (Y>-0.1) & (Y<0.1)); % position du deuxième chauffage (ht1)

k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A 

spy(G)                                % on dessine la grille
heat = G(H);                          % indices correspondant a H
heat1 =G(H1);
 
A = delsq(G);                         % matrice du -Laplacien sans 1/h^2 a l'interieur
door = [];                            % indices correspondant a la porte
window = [];                          % indices correspondant a la fenetre
for i = 2:n-1                        
  for j = 2:n-1
    ind = G(i,j);
    if ind ~= 0 
      if G(i,j+1) == 0                % la porte se situe sur le mur de droite
        if (Y(i,j)>0.3 && Y(i,j)<0.7) 
          door = [door ind];          % indices correspondant a la porte 
        else
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la porte le mur est isolant
        end
      end
      if G(i,j-1) == 0                % mur de gauche isolant: Neumann
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i-1,j) == 0                % mur en bas isolant: Neumann
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0                % la fenetre se situe sur le mur en haut
        if (X(i,j)>-0.5 && X(i,j)<-0.2) 
          window = [window G(i,j)];   % indices correspodant a la fenetre 
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
b(heat1) = -ht1;
b(heat) = -ht;                        % prise en compte du chauffage
b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*ot;                % CL Dirichet: fenetre
u = A\b;                              % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice 
                                      % correspondant a la grille
mesh(X,Y,U);                          % on dessine la solution sur la grille
axis('ij');
