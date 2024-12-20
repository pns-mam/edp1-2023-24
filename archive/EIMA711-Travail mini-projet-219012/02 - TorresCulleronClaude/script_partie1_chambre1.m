% PARTIE 1 CHAMBRE 1
clear all, close all
%Test �t� avec 20 degr�s pour la porte, la fen�tre et pas de chauffage
ot = 20;
dt = 20;
ht = 0;
n = 30;
titre = 'Temp�rature au primptemps, sans chauffage'
[U1, G1] = RoomTemperature1(ot,dt,ht,n,titre)

%Test en hiver avec 15 degr�s pour la porte, -10 pour la fen�tre 
% et pas de chauffage
ot = -10;
dt = 15;
ht = 0;
n = 30;
titre = 'Temp�rature en hiver, sans chauffage';
[U2, G2] = RoomTemperature1(ot,dt,ht,n,titre)

%Test en hiver avec 15 degr�s pour la porte, -10 pour la fen�tre
% et 300 de chauffage
ot = -10;
dt = 15;
ht = 300;
n = 30;
titre = 'Temp�rature en hiver, avec chauffage';
[U3, G3] = RoomTemperature1(ot,dt,ht,n,titre)


% Am�lioration de la position du chauffage
ot = -10;
dt = 15;
ht = 300;
n = 30;
titre = 'Temp�rature en hiver, avec le chauffage mieux plac�';
[U1, G1] = RoomTemperature3(ot,dt,ht,n,titre)