function U = TemperatureStationnaire2(ot,dt,ht,n)

%{
Cette fonction calculera la température ambiante dans une chambre
 Elle prend la température extérieure ot,
 la température de la porte dt et la température du chauffage ht et la
 nombre de points de grille n, puis calcule la température ambiante de la
 chambre
%}

x = linspace(-1,1,n);                            
[X,Y] = meshgrid(x,x); 
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
    
 
k = find(G);                        % les indices des G non nuls
G = zeros(size(G));                 % conversion logique - reel  
G(k) = (1:length(k))';              % suite d'indices de la matrice A  
A = delsq(G);                       % Laplacien a l'interieur du domaine G
door = []; 
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
u = A\b;                            % solution du systeme sous forme de vecteur
U = G;
U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
mesh(X,Y,U);                        % correspondant a la grille 
axis('ij'); 


disp("Température ambiante stationnaire "+mean(mean(u(G(G>0)))));