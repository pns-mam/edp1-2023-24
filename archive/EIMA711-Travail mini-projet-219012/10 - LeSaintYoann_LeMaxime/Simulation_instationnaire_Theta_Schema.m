n = 18;
ot= 35;       %Text
dt= 35;       %Tporte
ht = -2000;       %T chauffage

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
    if val~=0 %3
    
    if G(i+1,j)==0 % tout les murs face au 'nord'
        if (X(i,j)<=-0.15 && X(i,j)>=-0.95) %porte
        d=[d val];
        elseif (X(i,j)<=0.9 && X(i,j)>=0.2) % fenêtre 1
        w=[w val];
        else
        A(val,val)=A(val,val)-1;
        end
    end
    
    if (G(i,j+1)==0 || G(i,j-1)==0 || G(i-1,j)==0) %murs face à l'est sud et ouest
        A(val,val)=A(val,val)-1; % A(val,val)=A(val,val)-1
    end
    
    end
end
end

h = 2/(N-1);
A = -A/h^2; 
b = zeros(length(k),1);
heat = G(H);

b(w) = -(1/(h^2))*ot;              
b(heat) = -ht;                     
b(d) = -(1/(h^2))*dt;
u = A\b;                            
U = G;
U(G>0) = full(u(G(G>0)));

error = eps;
nu = 1; alpha = .5;
deltat = alpha/(2*nu)*h^2;
B = b*nu*deltat;
V_n=35*ones(size(A)(1));

i=1; 
T_MAX = 10^2;

M1=eye(size(A))-0.5*nu*deltat*A; % Mat à 'côté' de Un+1 : (I?(1??)??tA)
M2=eye(size(A))+0.5*nu*deltat*A; % Mat à 'côté' de Un : (I?(1??)??tA)
U_n = M1\(M2*V_n-B);

while (i<=T_MAX || i~=1) && mean(U_n)>=25
  i+=1;
  V_n = U_n;
  U_n = M1\(M2*V_n-B);
  
  P = G;
  U(G>0) = full(U_n(G(G>0))); 

  mesh(Y,X,U);
  xlabel('Y');
  ylabel('X');
  axis('ij');
  pause(0.05);
endwhile
