% PARTIE 1 CHAMBRE 2
clear all, close all
%Test �t� avec 20 degr�s pour la porte, la fen�tre et pas de chauffage
ot = 20;
dt = 20;
ht = 0;
n = 30;
titre = 'Temp�rature au primptemps, sans chauffage';
[U1, G1] = RoomTemperature2(ot,dt,ht,n,titre)

%Test en hiver avec 15 degr�s pour la porte, -10 pour la fen�tre 
% et pas de chauffage
ot = -10;
dt = 15;
ht = 0;
n = 30;
titre = 'Temp�rature en hiver, sans chauffage';
[U1, G1] = RoomTemperature2(ot,dt,ht,n,titre)

%Test en hiver avec 15 degr�s pour la porte, -10 pour la fen�tre
% et 300 de chauffage
ot = -10;
dt = 15;
ht = 300;
n = 30;
titre = 'Temp�rature en hiver, avec chauffage';
[U1, G1] = RoomTemperature2(ot,dt,ht,n,titre)

