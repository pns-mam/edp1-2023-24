L = 5;

% Crank-Nicolson
% Test de la precision spatiale
% On fixe cfl = 0.01, nt = 10 et on fait varier le pas d'espace
dx = L./[10 20 40 80];
err = [0.2381 0.074283 0.020812 5.3859e-03];
figure(1)
loglog(dx,err,'bx-',dx,dx.^2,'r*-'), grid on
title('Test de la precision spatiale')
legend('erreur','dx^2')
xlabel('dx')
ylabel('error')
