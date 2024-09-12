// Schema explicite pour l'equation d'advection
a=0.1;               // vitesse d'advection
Tmax=10;             // temps maximum
L=5;                 // longueur du domaine
deff('[v]=condinit(x)','v = bool2s((x>1.) & (x<1.5))');
N=201;               // Nombre points discretisation
dx=L/(N-1);          // pas de discretisation 
s=0.8;               // CFL <=1 assura la stabilite
dt=s*dx/a;           // calcul du pas de temps
x=linspace(0,L,N)'; 
u=condinit(x);
tps=dt:dt:Tmax; 
for t=tps            // boucle en temps
  uold=u; 
  u(2:N)=uold(2:N)-s*(uold(2:N)-uold(1:N-1)); 
  u(1)=u(N);
end 
t=tps($);
uexact=condinit(x-a*t);
plot(x,u,'rx-',x,uexact,'bo-');
xtitle('Comparaison entre la solution exacte et approchee au temps final','x','Solution');
hl=legend(['Sol. approchee';'Sol. exacte'],a=1);
