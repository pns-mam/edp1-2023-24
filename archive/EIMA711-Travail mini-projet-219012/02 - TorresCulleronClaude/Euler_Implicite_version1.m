function Euler_Implicite_version1(A,G,b,door,window,dt,ot,n,deltat, temp_init, nb_iter,affichage, timer)
  %% A,G,b,door,window,dt,wt sont les paramètres donnés par ... 
  %% n la discrétisation en espace
  %% deltat le pas de temps
  %% temp_init est la température initiale dans la chambre
  %% nb_iter le nombre de fois qu'on va parcourir la boucle while au maximum
  %% affichage est un booleen pour savoir si on affiche les graphiques
  %% timer est un booleen pour savoir si on chronomètre le temps pris a calculer
  h = 2/(n-1);
  alpha = 0.5; % condition CFL
  nu = 1;
  x = linspace(-1,1,n);                            
  [X,Y] = meshgrid(x,x);
  u = ones(length(b),1) * temp_init;
  u(door) = dt;
  u(window) = ot;

uold = zeros(length(b),1);
t = 0;
tStart = tic; 
while ( norm(uold-u) > 0.01 && t < nb_iter)
  t = t+1;
  uold = u;
  u = (eye(size(A)) - nu * deltat * A) \ (uold - nu * deltat * b);
  if(affichage)
    U = G;
    U(G>0) = full(u(G(G>0)));
    mesh(X,Y,U);
    title(t*deltat)
    axis('ij');
    pause(0.01);
  end
endwhile
tEnd = toc(tStart);
if(timer)
  disp('temps Euler Implicite : ' )
  disp(tEnd)
end  
