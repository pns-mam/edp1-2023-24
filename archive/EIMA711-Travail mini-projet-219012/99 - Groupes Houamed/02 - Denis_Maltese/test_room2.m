close all
figure('Name','Chambre 2:Stationnaire: porte � 20�')
room2(20,20,0,0,40);
figure('Name','Stationnaire: porte � 15�,-10� dehors sans chauffage')
room2(-10,15,0,0,40);
figure('Name','Stationnaire: porte � 15�,-10� dehors avec chauffage')
room2(-10,15,200,200,40);
disp("Partie 2")
disp("Cas Test 1:Solution stationnaire")

figure('Name','Solution stationnaire: retour de ski et chauffage ')
room2(-10,0,200,200,40);
disp("Cas Test 1:Euler explicite")
figure('Name','Euler explicite: retour de ski et chauffage ')
tic
room2eulexp(-10,0,200,200,40,1,0.5,0.1,-10);
toc
disp("Cas Test 1:Euler implicite")
figure('Name','Euler implicite: retour de ski et chauffage ')
tic
room2eulimp(-10,0,200,200,40,1,29,0.1,-10);
toc
disp("Cas Test 1:Crank-Nicholson")
figure('Name','Crank-Nicholson: retour de ski et chauffage ')
tic
room2cn(-10,0,200,200,40,1,29,0.1,-10);
toc
disp("Cas Test 2:Solution stationnaire")
figure('Name','Solution stationnaire: temp�rature id�ale �t� ')
room2(35,25,-100,-100,40);
disp("Cas Test 2:Euler explicite")
figure('Name','Euler explicite: temp�rature id�ale �t� ')
tic
room2eulexpId(35,25,-300,-300,40,1,0.5,0.1,35,18);
toc
disp("Cas Test 2:Euler implicite")
figure('Name','Euler implicite: temp�rature id�ale �t� ')
tic
room2eulimpId(35,25,-300,-300,40,1,29,1,35,18);
toc
disp("Cas Test 2:Crank-Nicholson")
figure('Name','Crank-Nicholson : temp�rature id�ale �t�')
tic
room2cnId(35,25,-300,-300,40,1,29,1,35,18);
toc

