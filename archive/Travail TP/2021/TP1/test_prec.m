
%% Test de la precision spatiale
L = 5;
dx = L./[10 20 40 80];
err = [0.1180 0.0348 0.0093 0.0024];
figure(1)
loglog(dx,err,'bx-',dx,dx.^2,'r*-'), grid on
title('Test de la precision spatiale')
legend('erreur','dx^2')
xlabel('dx')
ylabel('error')

%% Test de la precision en temps
dt = [0.0048 0.0019 9.6117e-04];
errt = [0.2424 0.1109 0.0580];
figure(2)
loglog(dt,errt,'bx-',dt,100*dt,'r*-'), grid on
title('Test de la precision en temps')
legend('erreur','dt')
xlabel('dt')
ylabel('error')
