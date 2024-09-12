////////////////////////////////////////////////////////////////
// Copyright G. Allaire, Juillet 2001, Janvier 2012
//
// Sur divers schemas de differences finies
// pour l'equation de transport:
// comparaison des schemas
//
////////////////////////////////////////////////////////////////
//
// equation de transport
// u,t + a u,x = 0
//
////////////////////////////////////////////////////////////////
//
lg = 1. ; // lg = longueur du domaine
dx = 0.01 ; // dx = pas d'espace
nx = lg/dx ; // nx = nombre de mailles
a = 1. ; // a = vitesse
cfl = 0.9 ; // cfl
dt = dx*cfl ; // dt = pas de temps
tfinal = 5. ; // tfinal = temps de simulation
nt = int(tfinal/dt) ; // nt = nombre de pas de temps effectues
//
// initialisation: condition initiale reguliere (sinus) 
//                                       ou singuliere (creneau)
//
x=zeros(1,nx) ;
u0=zeros(1,nx) ;
for i=1:nx
  x(i) = (i-1)*dx  ;
  if x(i) < 0.3
    u0(i) = -0.5 ;
  elseif x(i) > 0.7
     u0(i) = -0.5 ;
  else
     u0(i) = 0.5 ;
  end
//  u0(i) = sin(x(i)*2*%pi) ;
end
u=u0 ;
up=u0 ;
um=u0 ;
u1=u0 ;
up1=u0 ;
um1=u0 ;
u2=u0 ;
up2=u0 ;
um2=u0 ;
uexacte=u0 ;

tics=[4,10,4,10];
plotframe([0.,-1.1,lg,1.1],tics);
plot2d(x,u0,1,"000")
xtitle ('donnee initiale' ,' ',' ') ;


halt() ;
////////////////////////////////////////////////////////////////
// u = schema de Lax-Friedrichs
// u1 = schema de Lax-Wendroff
// u2 = schema decentre amont
////////////////////////////////////////////////////////////////
//
for n=1:nt
//
up = shiftp('+1',u) ;
um = shiftp('-1',u) ;

u=(up+um)/2 - a*dt/(2*dx)*(up-um) ;

up1 = shiftp('+1',u1) ;
um1 = shiftp('-1',u1) ;

u1=u1 - a*dt/(2*dx)*(up1-um1) + (0.5*dt*dt/(dx*dx))*a*a*( up1-2*u1+um1 ) ;

um2 = shiftp('-1',u2) ;

u2=u2 - a*dt/(dx)*(u2-um2) ;

if modulo(n,5) == 0
	for i=1:nx
	  uexacte(i) = -0.5 ;
	  arg(i) = x(i)-n*a*dt ;
              while arg(i) < 0 
                  arg(i) = arg(i) + lg ;
              end
	  if abs(arg(i)-0.5) <= 0.2
	     uexacte(i) = 0.5 ;
	  end
	end
//        for i=1:nx
//             uexacte(i) = sin((x(i)-a*n*dt)*2*%pi) ;
//        end
        clf()
        plotframe([0.,-1.1,lg,1.1],tics);
        plot2d(x,uexacte,[1,1],"100","solution exacte")
        plot2d(x,u,[2,2],"100","schema de Lax-Friedrichs")
        plot2d(x,u1,[3,3],"100","schema de Lax-Wendroff")
        plot2d(x,u2,[4,4],"100","schema decentre amont")
        xtitle('schemas explicites pour advection',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
	for i=1:nx
	  uexacte(i) = -0.5 ;
	  arg(i) = x(i)-nt*a*dt ;
              while arg(i) < 0 
                  arg(i) = arg(i) + lg ;
              end
	  if abs(arg(i)-0.5) <= 0.2
	     uexacte(i) = 0.5 ;
	  end
	end
//        for i=1:nx
//             uexacte(i) = sin((x(i)-a*nt*dt)*2*%pi) ;
//        end
        clf()
        plotframe([0.,-1.1,lg,1.1],tics);
        plot2d(x,uexacte,[1,1],"100","solution exacte")
        plot2d(x,u,[2,2],"100","schema de Lax-Friedrichs")
        plot2d(x,u1,[3,3],"100","schema de Lax-Wendroff")
        plot2d(x,u2,[4,4],"100","schema decentre amont")
        xtitle('schemas explicites pour advection, temps=5',' ',' ');

halt() ;
////////////////////////////////////////////////////////////////
// schema de Lax-Friedrichs
////////////////////////////////////////////////////////////////
u=u0 ;
u1=u0 ;
dt1=dt/2 ;
//
for n=1:nt
//
up = shiftp('+1',u) ;
um = shiftp('-1',u) ;

u=(up+um)/2 - a*dt/(2*dx)*(up-um) ;

up1 = shiftp('+1',u1) ;
um1 = shiftp('-1',u1) ;

u1=(up1+um1)/2 - a*dt1/(2*dx)*(up1-um1) ;

up1 = shiftp('+1',u1) ;
um1 = shiftp('-1',u1) ;

u1=(up1+um1)/2 - a*dt1/(2*dx)*(up1-um1) ;

if modulo(n,5) == 0
	for i=1:nx
	  uexacte(i) = -0.5 ;
	  arg(i) = x(i)-n*a*dt ;
              while arg(i) < 0 
                  arg(i) = arg(i) + lg ;
              end
	  if abs(arg(i)-0.5) <= 0.2
	     uexacte(i) = 0.5 ;
	  end
	end
//        for i=1:nx
//             uexacte(i) = sin((x(i)-a*n*dt)*2*%pi) ;
//        end
        clf()
        plotframe([0.,-1.1,lg,1.1],tics);
        plot2d(x,uexacte,[1,1],"100","solution exacte")
        plot2d(x,u,[2,2],"100","schema de Lax-Friedrichs cfl=0.9")
        plot2d(x,u1,[3,3],"100","schema de Lax-Friedrichs cfl=0.45")
        xtitle('influence de la cfl sur le schema de Lax-Friedrichs',' ',' ');
        xpause(100000) ;
end
//
end


// comparaison au temps final
	for i=1:nx
	  uexacte(i) = -0.5 ;
	  arg(i) = x(i)-n*a*dt ;
              while arg(i) < 0 
                  arg(i) = arg(i) + lg ;
              end
	  if abs(arg(i)-0.5) <= 0.2
	     uexacte(i) = 0.5 ;
	  end
	end
//        for i=1:nx
//             uexacte(i) = sin((x(i)-a*n*dt)*2*%pi) ;
//        end
        clf()
        plotframe([0.,-1.1,lg,1.1],tics);
        plot2d(x,uexacte,[1,1],"100","solution exacte")
        plot2d(x,u,[2,2],"100","schema de Lax-Friedrichs cfl=0.9")
        plot2d(x,u1,[3,3],"100","schema de Lax-Friedrichs cfl=0.45")
        xtitle('influence de la cfl sur le schema de Lax-Friedrichs, temps=5',' ',' ');

halt() ;
////////////////////////////////////////////////////////////////
// schema decentre amont
////////////////////////////////////////////////////////////////
u=u0 ;
u1=u0 ;
dt1=dt/2 ;
//
for n=1:nt
//
um = shiftp('-1',u) ;

u=u - a*dt/(dx)*(u-um) ;

um1 = shiftp('-1',u1) ;

u1=u1 - a*dt1/(dx)*(u1-um1) ;

um1 = shiftp('-1',u1) ;

u1=u1 - a*dt1/(dx)*(u1-um1) ;

if modulo(n,5) == 0
	for i=1:nx
	  uexacte(i) = -0.5 ;
	  arg(i) = x(i)-n*a*dt ;
              while arg(i) < 0 
                  arg(i) = arg(i) + lg ;
              end
	  if abs(arg(i)-0.5) <= 0.2
	     uexacte(i) = 0.5 ;
	  end
	end
//        for i=1:nx
//             uexacte(i) = sin((x(i)-a*n*dt)*2*%pi) ;
//        end
        clf()
        plotframe([0.,-1.1,lg,1.1],tics);
        plot2d(x,uexacte,[1,1],"100","solution exacte")
        plot2d(x,u,[2,2],"100","schema decentre amont cfl=0.9")
        plot2d(x,u1,[3,3],"100","schema decentre amont cfl=0.45")
        xtitle('influence de la cfl sur le schema decentre amont',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
	for i=1:nx
	  uexacte(i) = -0.5 ;
	  arg(i) = x(i)-n*a*dt ;
              while arg(i) < 0 
                  arg(i) = arg(i) + lg ;
              end
	  if abs(arg(i)-0.5) <= 0.2
	     uexacte(i) = 0.5 ;
	  end
	end
//        for i=1:nx
//             uexacte(i) = sin((x(i)-a*n*dt)*2*%pi) ;
//        end
        clf()
        plotframe([0.,-1.1,lg,1.1],tics);
        plot2d(x,uexacte,[1,1],"100","solution exacte")
        plot2d(x,u,[2,2],"100","schema decentre amont cfl=0.9")
        plot2d(x,u1,[3,3],"100","schema decentre amont cfl=0.45")
        xtitle('influence de la cfl sur le schema decentre amont, temps=5',' ',' ');



