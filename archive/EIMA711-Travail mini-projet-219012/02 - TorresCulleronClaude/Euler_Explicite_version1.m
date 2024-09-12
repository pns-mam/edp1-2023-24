function Euler_Explicite_version1(A,G,b,door,window,dt,ot,n, temp_init, nb_iter,affichage, timer)
  % A,G,b,door,window,dt,wt sont les paramètres donnés par ... 
  % temp_init est la température initiale dans la chambre
  % nb_iter le nombre de fois qu'on va parcourir la boucle for
  % affichage est un booleen pour savoir si on affiche les graphiques
  % timer est un booleen pour savoir si on chronomètre le temps pris a calculer
  h = 2/(n-1);
  alpha = 0.5; % condition CFL
  nu = 1;
  deltat = alpha / (2 * nu) * h^2;
  x = linspace(-1,1,n);                            
  [X,Y] = meshgrid(x,x);
  u = ones(length(b),1) * temp_init;
  u(door) = dt;
  u(window) = ot;
tStart = tic; 
for l = 1:nb_iter
  uold = u;
  u = (eye(size(A)) + nu * deltat * A) * uold - nu * deltat * b;
  if(affichage)
    U = G;
    U(G>0) = full(u(G(G>0)));
    mesh(X,Y,U);
    title(l*deltat)
    axis('ij');
    pause(0.07);
  end
end
tEnd = toc(tStart);
if(timer)
  disp('temps EE' )
  disp(tEnd)
end  

