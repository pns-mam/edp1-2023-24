%Calcul du nombre d'etapes necessaires avant d'atteindre la temperature
%ideale en hiver
num_chambre = 1;
ot = -10;
dt = 20;
ht = 400;
n = 30;
[A,G,b,door,window] = Param_partie2(num_chambre, ot,dt,ht,n);
temp_init = 10;
nb_iter = 100;



deltat = 0.6;
affichage= true;
timer = false;
figure
t = Euler_Implicite_version2(num_chambre,A,G,b,door,window,dt,ot,n,deltat, temp_init, nb_iter,affichage, timer)
disp('t')
disp(t)



%Calcul du nombre d'etapes necessaires avant d'atteindre la temperature
%ideale en ete

num_chambre = 2;
ot = 38;
dt = 20;
ht = -300;
n = 30;
[A,G,b,door,window] = Param_partie2(num_chambre, ot,dt,ht,n);
temp_init = 30;
nb_iter = 100;

deltat = 0.4;
affichage= true;
timer = false;
figure
t = Euler_Implicite_version2(num_chambre,A,G,b,door,window,dt,ot,n,deltat, temp_init, nb_iter,affichage, timer)
disp('t')
disp(t)