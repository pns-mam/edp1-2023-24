function t = Euler_Implicite_version2(num_chambre,A,G,b,door,window,dt,ot,n,deltat, temp_init, nb_iter,affichage, timer)
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
r2 = 101;
r1 = 200; 
tStart = tic; 
while ( norm(uold-u) > 0.01 && t < nb_iter &&  (r2<19 || r2>21) && r1 > 100 )
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
  if(num_chambre == 2)
     r2 = mean(U(U>0));
  end
  if(num_chambre == 1)
    r1 = size(U(U<19),1) - size(U(U==0),1);
  end
end
tEnd = toc(tStart);
if(timer)
  disp('temps Euler Implicite : ' )
  disp(tEnd)
end  
