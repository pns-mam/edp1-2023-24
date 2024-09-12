clear all, close all
%% Pour cette partie plusieurs test sont disponibles
%% TEST = 1 : Euler Explicite avec affichage des graphiques
%% TEST = 2 : Temps pris par Euler Explicite
%% TEST = 3 : Euler Implicite avec affichage des graphiques
%% TEST = 4 : Temps pris par Euler Implicite

TEST = 3;

num_chambre = 2; % 1 ou 2 suivant la chambre choisie
ot = -10; %température de la fenetre
dt = 15; %température de la porte
ht = 300; %température du chauffage
n = 30; %discrétisation

[A,G,b,door,window] = Param_partie2(num_chambre, ot,dt,ht,n);
temp_init = 10; %température initiale de la chambre
nb_iter = 100; %nombre d'itérations maximum réalisées dans les boucles

if(TEST == 1)
  affichage = true;
  timer = false;
  figure
  Euler_Explicite_version1(A,G,b,door,window,dt,ot, n,temp_init, nb_iter,affichage, timer);
end

if(TEST == 2)
  affichage = false;
  timer = true;
  Euler_Explicite_version1(A,G,b,door,window,dt,ot, n,temp_init, nb_iter,affichage, timer);
end

if(TEST == 3)
  deltat = 0.4;
  affichage = true;
  timer = false;
  figure
  Euler_Implicite_version1(A,G,b,door,window,dt,ot,n,deltat, temp_init, nb_iter,affichage, timer);
end


if(TEST == 4)
  deltat = 0.5 / (2 * 1) * (2/(n-1))^2;
  disp(deltat)
  affichage = false;
  timer = true;
  Euler_Implicite_version1(A,G,b,door,window,dt,ot,n,deltat, temp_init, nb_iter,affichage, timer);
end
