clear all;
close all;
%
%
%
% ATTENTION: la fonction peut-�tre tr�s longue � s'ex�cuter
%Tester avec un nombre_discretisations petit 
%et un n qui est la discretisation de d�part petit
%
function U = vitesse(ot,dt,ht,n,nombre_discretisations)

no=zeros(nombre_discretisations,1); %vecteurs qui va contenir les diff�rents pas de discr�tisation de l'espace

%initialisation des vecteurs qui vont contenir les temps d'ex�cutions des m�thodes
implicite=zeros(1,nombre_discretisations); 
explicite=zeros(1,nombre_discretisations);
  for t=1:nombre_discretisations %on compare pour plusieurs discr�tisations diff�rentes de l'espace
  x = linspace(-1,1,n);                            
  [X,Y] = meshgrid(x,x);                             % grille bi-dimensionnelle 
  G = (((X>-0.5) & (X<0.5) & (Y>-0.5) & (Y<0.5)) |...
      ((X>-1) & (X<-0.5) & (Y>-0.2) & (Y<0.2))|...
      ((X>0.5) & (X<1) & (Y>-0.2) & (Y<0.2)));
 
  H = ((X>-0.2) & (X<0.2) & (Y>-0.2) & (Y<0.2));      % chauffage au centre
  
  k = find(G);                        % les indices des G non nuls
  G = zeros(size(G));                 % conversion logique - reel  
  G(k) = (1:length(k))';              % suite d'indices de la matrice A  
  A = delsq(G);% Laplacien a l'interieur du domaine G
  door = [];                          % indices correspondant a la porte
  for i = 2:n-1                       % ajout de CL Neumann
    for j = 2:n-1                     % pour les murs isolants
      ind = G(i,j);
      if ind~=0
        if G(i,j-1) == 0              % CL Neumann sur le mur gauche
          if (Y(i,j)>-0.2 && Y(i,j)<0.2) 
            door = [door ind];
          else
          A(ind,ind) = A(ind,ind)-1;
          end
        end
        if G(i,j+1) == 0 && j<n-1         % CL Neumann sur le mur droit
          A(ind,ind) = A(ind,ind)-1;
        end
        if G(i+1,j) == 0      % On garde la CL Dirichlet pour la fenetre
          A(ind,ind) = A(ind,ind)-1;  % CL Neumann sur le mur du bas 
        end
        if G(i-1,j) == 0
            A(ind,ind) = A(ind,ind)-1; % CL Neumann sur le mur du haut
        end 
      end
    end
  end
  h = 2/(n-1);
  A = -A/h^2;                         % division par h^2
  b = zeros(length(k),1);
  window = G(G(:,end-1)>0,end-1);     % indices correspondant a la fenetre

  heat = G(H);                        % indices correspondant au radiateur faire attention qu'il se trouve bien dans la pi�ce
  b(window) = -1/h^2*ot;              % prise en compte des CL Dirichlet
  b(heat) = -ht;                      % et du chauffage
  b(door) = -1/h^2*dt;
  I = eye(size(A));

  alpha = 0.5; nu = 1;
  deltat = alpha/(2*nu)*(h^2);
  u = ot*ones(length(k),1); %initialisation de la temp�rature de la pi�ce 
  temps_total=10;    %on garde le meme temps pour chaque discr�tisation
  iterations=round(temps_total/deltat);     
  for methode=1:2        %on calcule pour les 2 m�thodes
    if methode==1        %euler explicite
      tic
        for l = 1:iterations
          uold = u;
          %u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b);
          u = (eye(size(A))+nu*deltat*A)*uold - nu*deltat*b;
          U = G;
          U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
          %if mod(l,100) == 0
            %mesh(X,Y,U);                        % correspondant a la grille 
            %axis('ij');                         % on dessine la solution sur la grille
            %pause(0.1);
          %end
        end                    % on dessine la solution sur la grille
       explicite(1,t)=toc;
       
       else
        tic
        for l = 1:iterations
          uold = u;
          u = (eye(size(A))-nu*deltat*A)\(uold - nu*deltat*b);
          U = G;
          U(G>0) = full(u(G(G>0)));           % on met la solution dans une matrice
          %if mod(l,100) == 0
            %mesh(X,Y,U);                        % correspondant a la grille 
            %axis('ij');                         % on dessine la solution sur la grille
            %pause(0.1);
          %end
        end                    
       implicite(1,t)=toc;
     
    end
  disp(explicite);
  disp(implicite);
end
disp(mean(mean(U((((X>-0.5) & (X<0.5) & (Y>-0.5) & (Y<0.5)) |...    % Pour avoir la temperature Moyenne uniquement sur notre chambre on lui redonne les dimensions de celle ci
    ((X>-1) & (X<-0.5) & (Y>-0.2) & (Y<0.2))|...                      %pour ne pas prendre en compte tout les 0 de la grille
    ((X>0.5) & (X<1) & (Y>-0.2) & (Y<0.2)))))));
no(t,1)=n;
n=n+20; %on passe a un pas de discretisation plus petit
end
%on affiche maintenant les graphes de la vitesse d'ex�cution
subplot(2,1,1);
scatter(no,explicite);
title('euler explicite');
xlabel('taille de la grille');
ylabel('temps en seconde');
subplot(2,1,2);
scatter(no,implicite);
title('euler implicite');
xlabel('taille de la grille');
ylabel('temps en seconde');


end

%exemple: n=20 puis n=40
v=vitesse(-10,8,170,20,2);
