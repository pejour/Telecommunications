clear
close all

Fe = 24000;  % fréquence d'échantillonnage
Rb = 3000;   % débit binaire 
Tb = 1/Rb;
Nb = 1000;
Te = 1/Fe;
fc = 4000;  % (BW)
N = 101;
bits = randi([0,1],1,Nb);  % Génération de Nb bits

% on utilisera les mêmes valeurs de Ts et Ns au cours du TP

Ts = Tb;     % durée entre symboles
Ns = Fe*Ts;  % nb d'échantillons utilisés par symbole


% réponse impulsionnelle d'un filtre passe-bas

hc = (2*fc/Fe)*sinc(2*fc*(-(N-1)*Te/2 :Te: (N-1)*Te/2));
Hc = abs(fft(hc));
Frequences = 1/Te*(1:length(Hc));


%% Filtres d'émission


% 1ère chaine de transmission

x1 = bits*2 - 1;

% réponse impulsionnelle
h1 = ones(1,Ns);
somme1 = kron(x1,[1 zeros(1,Ns-1)]);

NRZ1 = filter(h1,1,somme1);
figure(1);
plot(NRZ1);
title('Réponse en sortie du filtre émission sur la 1e chaine');


% 2ème chaine de transmission

alpha = 0.5;
x2 = bits*2 - 1;

% réponse impulsionnelle
h2 = rcosdesign(alpha,8,Ns);
somme2 = kron(x2,[1 zeros(1,Ns-1)]);

NRZ2 = filter(h2,1,[somme2 zeros(1,8*Ns)]);
figure(2);
plot(NRZ2);
title('Réponse en sortie du filtre émission sur la 2e chaine');


%% Filtres de réception

% 1ère chaine de transmission

x1 = bits*2 - 1;

% réponse impulsionnelle
hr1 = h1;

yc1 = filter(hc,1, [NRZ1 zeros(1,(N-1)/2)]);
zr1 = filter(hr1,1,yc1((N-1)/2+1:end));     % zeros(1,(N-1)/2) => yc1((N-1)/2+1:end) pour centrer le diagramme de l'oeil
figure(3);
plot(zr1);
title('Réponse en sortie du filtre de réception sur la 1e chaine');


% 2ème chaine de transmission

alpha = 0.5;
x2 = bits*2 - 1;

% réponse impulsionnelle
hr2 = h2;

yc2 = filter(hc,1,[NRZ2 zeros(1,(N-1)/2)]);
zr2 = filter(hr2,1,yc2((N-1)/2+1:end));
figure(4);
plot(zr2);
title('Réponse en sortie du filtre de réception sur la 2e chaine');



%% Echantillonnage

% Calcul de g1 et g2, fonctions de convolution
g1 = conv(h1,hr1);
g2 = conv(h2,hr2);

figure(5);
plot(g1);
title('Réponse impulsionnelle globale de la 1e chaine de transmission');

figure(6);
plot(g2);
title('Réponse impulsionnelle globale de la 2e chaine de transmission');



figure(7);
plot(Frequences,abs(fftshift(fft(g1,N))),Frequences,100*fftshift(Hc));
title("|Hc(f)|(orange) et |H(f)Hr(f)|(bleu) pour la chaine 1");
xlabel('Fréquences(Hz)');
ylabel('fft');

figure(8);
plot(Frequences,abs(fftshift(fft(g2,N))),Frequences,10*fftshift(Hc));
title("|Hc(f)|(orange) et |H(f)Hr(f)|(bleu) pour la chaine 2");
xlabel('Fréquences(Hz)');
ylabel('fft');

% pour les figures 7 et 8 on a centré chacune des courbes obtenues avec
% "fftshift" et on a ajusté la taille des portes


t0_1 = Ts;
n0_1 = t0_1/Te;
t0_2 = Te;
n0_2 = 1;   % = t0_2/Te


% tracés des diagrammes de l'oeil

% (une ligne apparait en plus, elle correspond au premier bit émis)

figure(9);
plot(reshape(zr1,Ns,length(zr1)/Ns));
title("Diagramme de l'oeil en sortie du filtre de la 1e chaine de transmission");

figure(10);
plot(reshape(zr2,Ns,length(zr2)/Ns));
title("Diagramme de l'oeil en sortie du filtre de la 2e chaine de transmission");


% signaux échantillonnés
z1_ech = zr1(1,(n0_1:Ns:end));
z2_ech = zr2(1,(Ns*8+n0_2:Ns:end));


%% Decision

d1 = sign(z1_ech);
d2 = sign(z2_ech);


%% Demapping

z_res_1 = (d1 + 1)/2;
z_res_2 = (d2 + 1)/2;


% calcul des taux d'erreur binaire
nb_bits_errones_1 = length(find(z_res_1-bits ~= 0));
nb_bits_errones_2 = length(find(z_res_2-bits ~= 0));
nb_bits_totaux = length(bits);

tx_erreur_binaire_1 = nb_bits_errones_1/nb_bits_totaux;
tx_erreur_binaire_2 = nb_bits_errones_2/nb_bits_totaux;



