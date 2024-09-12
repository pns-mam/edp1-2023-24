clf;
clear;
function Un=heat(xspan, tspan, nu, u0, f)
    Un = u0(1,2:length(u0)-1)';
    N=length(xspan)-2;
    A=2*diag(ones(N,1),0)-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1);
    I=diag(ones(N,1),0);
    sigma=nu*tspan(2)/(xspan(2)**2);
    for n = 1 : length(tspan)-1
        Fn=tspan(2)*f(xspan(1,2:N+1),tspan(n));
        Fn1=tspan(2)*f(xspan(1,2:N+1),tspan(n+1));
        B=(I-sigma*A/2)*Un+Fn'/2+Fn1'/2;
        Un1 = (I+sigma*A/2)\B;
        Un=Un1;
    end
    Un=[0;Un;0];
endfunction

function [y]=_f(xspan,t)
    y = -sin(xspan)*sin(t)+sin(xspan)*cos(t);
endfunction

function [y]=_u0(xspan)
    y=sin(xspan);
endfunction
function [y]=_uexact(xspan,t)
    y=sin(xspan)*cos(t);
endfunction

L=%pi;
N=50;
dx=L/(N-1);
s=0.5;
nu=1
dt=s*dx*dx/nu;
tmax=1;
//N=50;
x=linspace(0,L,N);
t=0:dt:tmax
//u0=_u0(dx:dx:%pi-dx);
u0=_u0(x);
u=heat(x,t,nu,u0,_f);
uexact=_uexact(x,1);
plot(x,u,'rx-',x,uexact,'bo-');
xtitle('Comparaison entre la solution exacte et approchee', 'x', 'Solution');
legend('Sol. approchee', 'Sol. exacte');