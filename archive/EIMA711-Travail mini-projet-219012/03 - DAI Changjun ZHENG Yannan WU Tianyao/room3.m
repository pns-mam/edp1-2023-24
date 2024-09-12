
ot=-10;
dt=15;
ht=300;
n=40;

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
G = ((X<0)&X.^2+(Y-0.2).^2<0.64 | ...
    (X>0 & X.^2+(Y+0.2).^2<0.64&Y>0)|...
    (X>0&Y<0&X<0.6&Y>-0.6));
                  
H = ((X>-0.3) & (X<0.3) & (Y>-0.2) & (Y<0.4));      % chauffage au centre
% H = ((X-0.1).^2+(Y+0.1).^2<0.1);     % chauffage entre deux
%H = ((X).^2+(Y).^2<0.1);  
% pres de la mur

k = find(G);                          % les indices des G non nuls
G = zeros(size(G));                   % conversion logique - reel
G(k) = (1:length(k))';                % suite d'indices de la matrice A
%spy(G)                                % on dessine la grille
heat = G(H);                          % indices correspondant a H

A = delsq(G);                         % matrice du -Laplacien sans 1/h^2 a l'interieur
door = [];                            % indices correspondant a la porte
window = [];                          % indices correspondant a la fenetre
for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
        if (X(i,j)>=0&&Y(i,j)>0)
            if(G(i-1,j) == 0 ||G(i+1,j) == 0 || G(i,j+1) == 0)
                     A(ind,ind) = A(ind,ind)-1; 
            end
        end
        if (X(i,j)>0&&Y(i,j)<0)
            if G(i-1,j)==0
                door=[door ind];
            end
            if G(i,j+1)==0
                A(ind,ind) = A(ind,ind)-1;
            end
        end
        if (X(i,j)<0)
            if(G(i-1,j) == 0 ||G(i+1,j) == 0 || G(i,j-1) == 0)
                    window = [window G(i,j)];
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
u = A\b;                        % solution du systeme sous forme de vecteur                              
U = G;
U(G>0) = full(u(G(G>0)));             % on met la solution dans une matrice 
                                      % correspondant a la grille
mesh(X,Y,U);                          % on dessine la solution sur la grille
view(50,50);
axis('ij');