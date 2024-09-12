function U = poisson_instat1_imp(ht, dt, wt, n) 
% Simulation instationaire du probleme de Poisson
% Version avec schema d'euler implicite, chambre 1
%
% ht temperature du chauffage
% dt temperature de la porte d'entree
% wt temperature de la fenetre
% n nombre de points de discretisation
U_exact = poisson_stat(ht, dt, wt, n); %solution exacte pour vérifier la convergence

x = linspace(-1,1,n);
[X,Y] = meshgrid(x,x);                  % grille bi-dimensionnelle
G = ((X>-0.9) & (X<0) & (Y>-0.7) & (Y<0.8)) | ...  % geometrie de la chambre 1
   ((X>0) & (X<0.7) & (Y>-0.4) & (Y<0.8)) | ... 
   ((X>0.7) & (X<1) & (Y>0.4) & (Y<0.8));

H = ((X>-0.8) & (X<-0.5) & (Y>0.7) & (Y<0.8));       % position du chauffage
%H = ((X>-0.3) & (X<0.2) & (Y>0.1) & (Y<0.3));       % autre option de terme de source

k = find(G);                          % les indices des G non nuls
G = zeros(size(G));                   % conversion logique - reel
G(k) = (1:length(k))';                % suite d'indices de la matrice A
spy(G)                                % on dessine la grille
heat = G(H);                          % indices correspondant a H

A = delsq(G);                         % matrice du Laplacien dans grille
door = [];
window = [];
% pour les commentaires de cette initialisation,
% cf la fonction poisson_stat1
for i = 2:n-1                        
  for j = 2:n-1
    ind = G(i,j);
    if ind ~= 0 
      if G(i-1,j) == 0                
        if (X(i,j)>0.1 && X(i,j)< 0.6)
          window = [window G(i,j)];   
        else
          A(ind,ind) = A(ind,ind)-1;  
        endif
      endif
      if G(i+1,j) == 0                
        A(ind,ind) = A(ind,ind)-1;
      endif
      if G(i,j-1) == 0                
        A(ind,ind) = A(ind,ind)-1;
      endif
      if G(i,j+1) == 0     
        if (Y(i,j)>0.5 && Y(i,j)<0.7) 
          door = [door ind];          
        else
          A(ind,ind) = A(ind,ind)-1;
        endif
      endif
    endif
  endfor
endfor

h = 2/(n-1);                          % pas du maillage
A = -A/h^2;                           % division par h^2
b = zeros(length(k),1);               % initialisation du second membre
b(heat) = -ht;                        % prise en compte du chauffage
b(door) = -1/h^2*dt;                  % CL Dirichlet: porte
b(window) = -1/h^2*wt;                % CL Dirichet: fenetre

U = G;
nb_step = 0;
alpha = 0.5; 
nu =1;
deltat = alpha/(2*nu)*h^2;
u = 40*ones(length(k),1);     % température initiale uniforme
tic;
while nb_step < 300 %critère d'arrêt à étapes arbitraire
%while abs(sum((U-U_exact)(:))) > 1 %critère d'arrêt à précision arbitraire             
    nb_step++;
    uold = u;
    u = (eye(size(A))-nu*deltat*A)\(uold - deltat*b); % implicite
    U = G;
    U(G>0) = full(u(G(G>0)));
    mesh(X,Y,U);
    axis('ij');
    title(sprintf("Simulation par schema d'euler implicite, t=%d", nb_step));
    pause(0.01);
endwhile
time = toc