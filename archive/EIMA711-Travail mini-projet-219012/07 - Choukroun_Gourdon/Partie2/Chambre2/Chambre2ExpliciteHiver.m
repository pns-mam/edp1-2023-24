ot = -10;
dt = 15;
ht = 430;
n = 30;

%{
Code permettant la simulation de chauffage d'une pièce en hiver avec une température extérieure à -10°C
une température à la porte de 15°C, et le chauffage apportant de la chaleur
à une puissance de 430
On utilise ici le schéma explicite.
%}

x = linspace(-1,1,n);   % x : vecteur de -1 à 1 avec n valeurs équidistantes
% X : matrice avec les valeurs de x en ligne donc toutes ses lignes sont égales à x
% Y : matrice avec les valeurs de x en colonnes donc toutes ses colonnes sont égales à x
[X,Y] = meshgrid(x,x); % grille bi-dimensionnelle

G=(((X>-0.2) & (X<1) & (Y>-0.5) & (Y<1)) | ... %pièce principale
    ((X>0.2) & (X<0.7) & (Y>-0.8) & (Y<-0.5))  | ... %coin détente
    ((X>-0.5) & (X<-0.2) & (Y>0.1) & (Y<0.5)) | ...  %arche,couloir
    ((X>-0.8) & (X<-0.5) & (Y>0) & (Y<0.7)) |... %salle de bain
    ((X>-1) & (X<-0.8) & (Y>0.3) & (Y<0.5)) | ... %fenêtre de la salle de bain
    ((X>-0.1) & (X<0.1) & (Y>-0.7) & (Y<-0.5)) | ... %fenêtre de la pièce principale
    ((X>0.3) & (X<0.55) & (Y>-1) & (Y<-0.8))) %fenêtre du coin détente

%         CHAUFFAGES    

%¨1e positionnement des chauffages (2 positions  : un sous la
%fenêtre de la salle de bain l'autre sous la fenêtre du coin détente) avec
%des tailles 'faibles'
%H = (((X>-0.8) & (X<-0.7) & (Y>0.3) & (Y<0.5)) | ...
  %  ((X>0.3) & (X<0.5) & (Y>-0.8) & (Y<-0.7))) 
  
 %¨2e positionnement des chauffages (2 positions  : un sous la
%fenêtre de la salle de bain l'autre sous la fenêtre de la pièce
%principale) avec des tailles 'faibles'
%H = (((X>-0.8) & (X<-0.7) & (Y>0.3) & (Y<0.5)) | ...
   % ((X>-0.1) & (X<0.1) & (Y>-0.5) & (Y<-0.4)))
   
   %¨3e positionnement des chauffages (2 positions  : un sous la
%fenêtre de la salle de bain l'autre sous la fenêtre du coin détente)-
%mêmes que le 1e avec un chauffage de plus grande taille pour augmenter son
%efficacité => CHAUFFAGE OPTIMAL
H = (((X>-0.8) & (X<-0.6) & (Y>0.15) & (Y<0.58)) | ...
    ((X>0.2) & (X<0.59) & (Y>-0.8) & (Y<-0.58)))  

   %¨4e positionnement des chauffages (2 positions différentes : un sous la
%fenêtre de la salle de bain l'autre sous la fenêtre de la pièce
%principale)- mêmes que le 3e avec un chauffage de plus grande taille pour augmenter son
%efficacité
%H = (((X>-0.8) & (X<-0.6) & (Y>0.15) & (Y<0.58)) | ...
 %   ((X>-0.19) & (X<0.2) & (Y>-0.5) & (Y<-0.28)))

%5e position un grand chauffage au sol au centre de la piece
    %rectangulaire, axé vers la parte où il y a les 2 fenêtres
    %H = ((X>-0.1) & (X<0.6) & (Y>-0.5) & (Y<-0.2)) 

k = find(G);                        % on trouve les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel % G garde sa taille mais tous ses éléments deviennent nuls 
G(k) = (1:length(k))';              % suite d'indices de la matrice A % G reprend sa forme de carré central avec cette fois les valeurs du carré allant de 1 au nombre de valeurs, s'implémentant dans le sens des colonnes 
A = delsq(G);                       % matrice du Laplacien à l'intérieur du domaine G
door = [];                          % indices correspondant à la porte
windowSDB=[];
windowpiece=[];
windowcoin=[];

for i = 2:n-1                       % ajout de CL Neumann
  for j = 2:n-1                     % pour les murs isolants
    ind = G(i,j);
    if ind~=0
      if G(i,j-1) == 0    %murs de gauche          
        if (Y(i,j)>0.3 && Y(i,j)<0.5) %position de la fenêtre de la salle de bain
          windowSDB = [windowSDB ind];        % On garde la CL Dirichlet pour la fenêtre 
        else
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur de gauche
        end
      end
      if G(i,j+1) == 0              % CL Neumann sur le mur de droite
        A(ind,ind) = A(ind,ind)-1;
      end
      if G(i+1,j) == 0      % murs du bas
        if (X(i,j)>0 && X(i,j)<0.2) %position de la porte
          door = [door ind];        % On garde la CL Dirichlet pour la porte 
        else
          A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du bas
        end
      end
      if G(i-1,j) == 0      % murs du haut
        if (X(i,j)>-0.1 && X(i,j)<0.1) %position de la fenêtre de la pièce principale
          windowpiece = [windowpiece ind];        % On garde la CL Dirichlet pour la fenêtre
        else
            if (X(i,j)>0.3 && X(i,j)<0.55) %position de la fenêtre du coin détente
                windowcoin = [windowcoin ind];        % On garde la CL Dirichlet pour la fenêtre
            else
                A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
            end
         end
      end
    end
  end
end

h = 2/(n-1);
A = -A/h^2;                         % division par h^2
b = zeros(length(k),1);
heat = G(H);                        % indices correspondant au radiateur
b(windowSDB) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(windowpiece) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(windowcoin) = -1/h^2*ot;              % prise en compte des CL Dirichlet
b(heat) = -ht;                      % et du chauffage
b(door) = -1/h^2*dt;





[U,u] = Chambre2(ot,dt,0,n); %Condition initiale : on réalise l'approximation en stationnaire sans chauffage avec la même pièce
[U_final,u_final] = Chambre2(ot,dt,ht,n); % solution finale en stationnaire

% Paramètres de la discrétisation
L = 5 ;            % longueur du domaine
nx = 50 ;          % nombre de mailles
dx = L/(nx+1) ;    % pas d'espace
cfl = 0.1 ;        % cfl
nu = 1;            % coefficient de diffusion
dt2 = (dx^2/2*nu)*cfl ; % dt = pas de temps % car cfl = nu*dt/dx^2
nt = 10 ;          % nombre de pas de temps effectues

[m,n]=size(u);
I = eye(m);

tic
cpt=0
while (mean(mean(u(G(G>0)))) < 21)
    cpt=cpt+1;
    u_old = u;
    u=(I+nu*A*dt2)*u_old-b*dt2; %explicite
    U = G;
    U(G>0) = full(u(G(G>0)));       
    mesh(X,Y,U);                         
    axis('ij'); 
    pause(0.001)
end
toc



disp("Température méthode explicite "+mean(mean(u(G(G>0)))));
disp("Température en stationnaire "+mean(mean(u_final(G(G>0)))));
%disp("Compteur "+cpt);
