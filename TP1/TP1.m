Fe = 24000;
Rb = 6000;
Tb = 1/Rb;
Nb = 1000;

bits = randi([0,1],1,Nb);  % Génération de Nb bits

%% Modulateur 1

Ts1 = Tb;
Ns1 = Fe*Ts1;
x1 = bits*2 - 1;
h1 = ones(1,Ns1);
somme1 = kron(x1,[1 zeros(1,Ns1-1)]);

NRZ1 = filter(h1,1,somme1);
figure(1);
plot(NRZ1);

dsp1 = pwelch(NRZ1,[],[],[],Fe,'twosided');


% Calcul de la DSP théorique

figure(2);
Frequences1 = Fe/length(dsp1)*(1:length(dsp1));
semilogy(Frequences1, dsp1, Frequences1, (Ts1*(sinc(Frequences1*Ts1).^2) + Ts1*(sinc((Frequences1-Fe)*Ts1).^2)));


%% Modulateur 2

Ts2 = Tb*2;
Ns2 = Fe*Ts2;
x2 = 2*bi2de(reshape(bits,Nb/2,2))' - 3;

h2 = ones(1,Ns2);
somme2 = kron(x2,[1 zeros(1,Ns2-1)]);

NRZ2 = filter(h2,1,somme2);
figure(3);
plot(NRZ2);

dsp2 = pwelch(NRZ2,[],[],[],Fe,'twosided');

% Calcul de la DSP théorique (écart-type = 5)

figure(4);
Frequences2 = Fe/length(dsp2)*(1:length(dsp2));
semilogy(Frequences2, dsp2, Frequences2, 5*(Ts2*(sinc(Frequences2*Ts2).^2) + Ts2*(sinc((Frequences2-Fe)*Ts2).^2)));


%% Modulateur 3

Ts3 = Tb;
Ns3 = Fe*Ts3;
x3 = bits*2 - 1;
h3 = [-ones(1, Ns3/2) ones(1,Ns3/2)];

somme3 = kron(x3,[1 zeros(1,Ns3-1)]);

NRZ3 = filter(h3,1,somme3);
figure(5);
plot(NRZ3);

dsp3 = pwelch(NRZ3,[],[],[],Fe,'twosided');


% Calcul de la DSP théorique

figure(6);
Frequences3 = Fe/length(dsp3)*(1:length(dsp3));
dsp3th = Ts3*(sin(pi*Frequences3*Ts3/2).^4)./(pi*Frequences3*Ts3/2).^2;
DSP3 = dsp3th + flip(dsp3th);
semilogy(Frequences3, dsp3, Frequences3, DSP3);



%% Modulateur 4

Ts4 = Tb;
Ns4 = Fe*Ts4;
alpha = 0.5;
x4 = bits*2 - 1;
h4 = rcosdesign(alpha,8,Ns4);

somme4 = kron(x4,[1 zeros(1,Ns4-1)]);

NRZ4 = filter(h4,1,somme4);
figure(7);
plot(NRZ4);

dsp4 = pwelch(NRZ4,[],[],[],Fe,'twosided');


% sigma_a = 1 car symboles binaires à moyenne nulle

Frequences4 = Fe/length(dsp4)*(1:length(dsp4));

n1 = length(dsp4)*(1-alpha)/(2*Ts4*Fe);     % on passe par les indices où f = (1-alpha)/(2*Ts4*Fe) et
n2 = length(dsp4)*(1+alpha)/(2*Ts4*Fe);     % f = (1+alpha)/(2*Ts4*Fe)

dsp4_th = zeros(length(Frequences4));

dsp4_th(1:n1) = 1;
dsp4_th(n1:n2) = 1/2*(1+cos((pi*Ts4/alpha)*(Frequences4(n1:n2) - (1-alpha)/(2*Ts4))));

figure(8);
semilogy(Frequences4, dsp4/mean(dsp4(1:n1)), Frequences4, dsp4_th + flip(dsp4_th));



%% Afficher toutes les DSP pour les comparer

figure(9);
semilogy(Frequences1, dsp1/max(dsp1), Frequences2, dsp2/max(dsp2), Frequences3, dsp3/max(dsp3), Frequences4, dsp4/max(dsp4));



%% Réponses aux questions :


% 1) Quel est le classement des modulateurs bande de base étudiés par ordre d’efficacité spectrale croissante ?
% La courbe du modulateur 3 (en jaune) est celle où la largeur de la bande 
% est la plus grande, c'est donc le modulateur 3 qui est le moins efficace.
% Ensuite, il s'agit de la courbe bleue correspondant au modulateur 1
% (toujours sur le même critère de la largeur de la bande).
% Puis c'est la courbe du modulateur 2 (en orange), sa bande étant moins
% large que celle du modulateur 1.
% Enfin, la bande de la courbe du modulateur 4 (en violet) est la moins
% large, c'est donc le modulateur dont l'efficacité est la plus grande.



% 2) Quels sont les éléments d’un modulateur bande de base qui agissent sur l’efficacité spectrale obtenue pour la transmission et de quelle manière ?
% Pour agir sur l'efficacité il faut agir sur la largeur de la bande => il
% faut diminuer la puissance et la largeur de la bande pour améliorer
% l'efficacité.

