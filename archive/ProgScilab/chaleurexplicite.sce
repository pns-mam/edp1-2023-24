// Schema explicite pour l'equation de la chaleur
nu=5;                // coefficient de viscosite
Nbmax=500;           //nombre maximum de pas de temps
L=1;                 // longueur du domaine
deff('[v]=condinit(x)','v = bool2s((x>-1.) & (x<1.5))');
N=51;               // Nombre points discretisation
dx=L/(N-1);          // pas de discretisation 
s=0.5;               // CFL <=1/2 assure la stabilite
dt=s*dx*dx/nu;       // calcul du pas de temps
x=linspace(0,L,N)'; 
u=condinit(x); 
for n=1:Nbmax        // boucle en temps
  uold=u;
  u(1)=0; u(N)=0; 
  u(2:N-1)=(1-2*s)*uold(2:N-1)+s*(uold(3:N)+uold(1:N-2));
  if modulo(n,5) == 0
    clf()
    plot(x,u,'rx-'); 
  end
end 