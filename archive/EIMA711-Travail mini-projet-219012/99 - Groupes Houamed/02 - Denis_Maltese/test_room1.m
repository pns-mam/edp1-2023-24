close all

disp('stationnaire');
disp('Test : été');
figure('Name','stationnaire ete')
Room1Temperature(20,20,0,40);
disp('Test : hiver; sans chauffage');
figure('Name','stationnaire hiver sans chauffage')
Room1Temperature(-10,15,0,40);
disp('Test : hiver; sans chauffage');
figure('Name','stationnaire hiver avec chauffage')
Room1Temperature(-10,15,600,40);

disp('Test : instationnaire hiver')
disp('Euler explicite; chambre froide');
ti = 0; % température initiale
figure('Name','explicite chambre froide')
tic
Instationnaire1(-10,0,600,40,ti,"e",0.5,false);
toc
disp('Euler implicite; chambre froide');
figure('Name','implicite chambre froide')
tic
Instationnaire1(-10,0,600,40,ti,"i",25,false);
toc
disp('Crank-Nicolson; chambre froide');
figure('Name','crank-nicolson chambre froide')
tic
Instationnaire1(-10,0,600,40,ti,"cn",25,false);
toc
disp("Comparaison avec la simulation stationnaire");
figure('Name','stationnaire chambre froide')
Room1Temperature(-10,0,600,40);

disp('Test : instationnaire été')
disp('Euler explicite; chambre chaude');
ti = 30; % température initiale
figure('Name','explicite chambre chaude')
tic
Instationnaire1(30,30,-500,40,ti,"e",0.5,true);
toc
disp('Euler implicite; chambre chaude');
figure('Name','implicite chambre chaude')
tic
Instationnaire1(30,30,-500,40,ti,"i",25,true);
toc
disp('Crank-Nicolson; chambre chaude');
figure('Name','crank-nicolson chambre chaude')
tic
Instationnaire1(30,30,-500,40,ti,"cn",25,true);
toc
