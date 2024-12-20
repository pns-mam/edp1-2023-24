n=50;
ot=-10;       %Text
dt=20;       %Tporte
ht=1500;       %T chauffage

  N=2*n+1;
  X=linspace(-1,1,N);
  [X,Y]=meshgrid(X,X);
  
  %forme de la chambre en L
  
  G=(X>-1 & X<1 & Y>-1 & Y<=0) | (Y>=0 & Y<1 & X>-1 & X<=0); %forme de la chambre en L (upside-down)
  
  H = ((X>-1) & (X<-0.7) & (Y>-1) & (Y<-0.7)); %coordonnées du radiateur
  
  k = find(G);                        % les indices des G non nuls
  G = zeros(size(G));                 % conversion logique - reel  
  G(k) = (1:length(k))';              % suite d'indices de la matrice A   
  A = delsq(G);                       % Laplacien a l'interieur du domaine G
  d = [];
  w = [];
  
  for i=1:N-1
    for j=1:N-1
      val=G(i,j);
      if val~=0
        
        if G(i+1,j)==0 % tout les murs face au 'nord'
          if (X(i,j)<=-0.15 && X(i,j)>=-0.95) %porte
            d=[d val];
          elseif (X(i,j)<=0.9 && X(i,j)>=0.2) % fenêtre 1 
            w=[w val];
          else
            A(val,val)=A(val,val)-1;
          end
        end
        if G(i,j+1)==0 
          if ( Y(i,j)<-0.3 && Y(i,j)>-0.7 )
            w=[w val];
          else A(val,val)=A(val,val)-1;
          end
        end
        if (G(i,j-1)==0 || G(i-1,j)==0) %murs face à l'est sud et ouest
          A(val,val)=A(val,val)-1; % A(val,val)=A(val,val)-1
        end
        
      end
    end
  end
 
  h = 2/(N-1);
  A = -A/h^2; 
  b = zeros(length(k),1);
  heat = G(H);
  
  %W= (Y==0 & X<=0.7 & X>=0.5); %coo de fenêtre 1
  %d= (Y==(1-2/(n-1)) & X<=-0.1 & X>=-0.5); %coo de porte
  
  b(w) = -(1/(h^2))*ot;              
  b(heat) = -ht;                     
  b(d) = -(1/(h^2))*dt;
  u = A\b;                            
  U = G;
  U(G>0) = full(u(G(G>0)));
  
  

  mesh(Y,X,U);
  xlabel('Y');
  ylabel('X');
  axis('ij');     

  moy=mean(mean(U))