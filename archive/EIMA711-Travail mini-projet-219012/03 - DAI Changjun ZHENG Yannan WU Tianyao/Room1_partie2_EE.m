ht = 300;
dt = 15;
ot = -10;
n = 41;  % nombre de points de discretisation

x = linspace(-1,1,n);
[X,Y] = meshgrid(x,x);                % grille bi-dimensionnelle
G = ((X>-1) & (X<=0) & (Y>-1) & (Y<0.5) | ...
    (X>=0 & X.^2+(Y+0.2).^2<0.6) | ...
    (X>-0.9) & (X<=0) & (Y>=0.5) & (Y<1));
H = ((X-0.3).^2+(Y+0.2).^2<0.2);                 % position du chauffage

k = find(G);                          % les indices des G non nuls
G = zeros(size(G));                   % conversion logique - reel
G(k) = (1:length(k))';                % suite d'indices de la matrice A
%spy(G)                                % on dessine la grille
heat = G(H);                          % indices correspondant a H

A = delsq(G);                         % matrice du -Laplacien sans 1/h^2 a l'interieur
door = [];                            % indices correspondant a la porte
window = [];                          % indices correspondant a la fenetre
for i = 2:n-1
  for j = 2:n-1
    ind = G(i,j);
    if ind ~= 0
      if G(i,j-1) == 0                % la porte se situe sur le mur a gauche
        if (Y(i,j)>0.6 && Y(i,j)<0.9)
          door = [door ind];          % indices correspondant a la porte
        else
          A(ind,ind) = A(ind,ind)-1;  % en dehors de la porte le mur est isolant
        end
      end
      if X(i,j)>0
        if(G(i-1,j) == 0 ||G(i+1,j) == 0 || G(i,j+1) == 0)
          window = [window G(i,j)];   % indices correspodant a la fenetre
        end
      else
        if G(i+1,j) == 0                % mur du bas isolant: Neumann
          A(ind,ind) = A(ind,ind)-1;
        end
        if G(i-1,j) == 0                % mur du haut isolant: Neumann
          A(ind,ind) = A(ind,ind)-1;
        end
        if G(i,j+1) == 0
          A(ind,ind) = A(ind,ind)-1;
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
b(window) = -1/h^2*ot;                % CL Dirichet: fenetre
u = zeros(length(k),1);

U = G;
U(G>0) = full(u(G(G>0))); 
time=0;
space=U(3:29,3:18);                %Zone principale

alpha = 0.5; nu =1;
deltat = alpha/(2*nu)*h^2;

while mean(space(:))<15           %Chauffer la zone principale ¨¤ une temp¨¦rature confortable
  time=time+1;
  u=u+deltat*nu*(A*u-b);
  U(G>0) = full(u(G(G>0))); 
  mesh(X,Y,U);
  title(['At time=',num2str(time)]);
  view(50,50);
  axis('ij');
  pause(0.01);
  space=U(3:29,3:18);  
end



