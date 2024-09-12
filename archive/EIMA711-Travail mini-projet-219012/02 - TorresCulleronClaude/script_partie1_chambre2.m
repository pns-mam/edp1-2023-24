% PARTIE 1 CHAMBRE 2
clear all, close all
%Test été avec 20 degrés pour la porte, la fenêtre et pas de chauffage
ot = 20;
dt = 20;
ht = 0;
n = 30;
titre = 'Température au primptemps, sans chauffage';
[U1, G1] = RoomTemperature2(ot,dt,ht,n,titre)

%Test en hiver avec 15 degrés pour la porte, -10 pour la fenêtre 
% et pas de chauffage
ot = -10;
dt = 15;
ht = 0;
n = 30;
titre = 'Température en hiver, sans chauffage';
[U1, G1] = RoomTemperature2(ot,dt,ht,n,titre)

%Test en hiver avec 15 degrés pour la porte, -10 pour la fenêtre
% et 300 de chauffage
ot = -10;
dt = 15;
ht = 300;
n = 30;
titre = 'Température en hiver, avec chauffage';
[U1, G1] = RoomTemperature2(ot,dt,ht,n,titre)

