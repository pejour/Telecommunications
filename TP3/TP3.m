clear
close all

Fe = 24000;  % fréquence d'échantillonnage
Rb = 6000;   % débit binaire 
Tb = 1/Rb;
Nb = 10000;
Te = 1/Fe;
N = 101;
bits = randi([0,1],1,Nb);  % Génération de Nb bits

% on réutilise la 1ère chaine de transmission du TP2

Ts = Tb;     % durée entre symboles
Ns = Fe*Ts;  % nb d'échantillons utilisés par symbole
M = 2;  % ordre de la modulation

x1 = bits*2 - 1;
t0_1 = Ts;
n0_1 = t0_1/Te;


%% 2. Chaine de référence


% réponse impulsionnelle
h1 = ones(1,Ns);
somme1 = kron(x1,[1 zeros(1,Ns-1)]);

NRZ1 = filter(h1,1,somme1);
figure(1);
plot(NRZ1);
title('Réponse en sortie du filtre émission sur la chaine de référence');


%% Introduction du bruit

TEB = zeros(1,9);
TEBth =  zeros(1,9);


for i=0:8
    Px = mean(abs(NRZ1).^2);        % puissance du signal              
    SNRdb = i;
    SNR = 10^(SNRdb/10);       % SNR = EB/n0
    sigma = sqrt(Px*Ns/(2*log2(M)*SNR));

    bruit = sigma*randn(1,length(NRZ1));    % expression du bruit

    NRZ1_bruite = NRZ1 + bruit;
    
    zr1 = filter(h1,1,NRZ1_bruite);
    
    
    % signaux échantillonnés
    z1_ech = zr1(1,(n0_1:Ns:end));
    
    % décision
    d1 = sign(z1_ech);
    z_res_1 = (d1 + 1)/2;
    
    % calcul des taux d'erreur binaire
    nb_bits_errones_1 = length(find(z_res_1-bits ~= 0));
    nb_bits_totaux = length(bits);

    tx_erreur_binaire_1 = nb_bits_errones_1/nb_bits_totaux;
    
    TEB(i+1) = tx_erreur_binaire_1;
    
    TEBth(i+1) = qfunc(sqrt(2*SNR));
end


% tracé TEB théorique et expérimental en fonction du SNRdb

figure(2);
semilogy([0:8],TEB,[0:8],TEBth);
title('Evolution du TEB et du TEB théorique en fonction de SNRdb');

% on choisit initialement Nb = 10000 pour avoir une précision du TEB à 10%
% (voir annexe)



%% 3. Première chaine


h2 = [ones(1, Ns/2) zeros(1,Ns/2)];

somme2 = kron(x1,[1 zeros(1,Ns-1)]);

NRZ2 = filter(h1,1,somme2);
figure(3);
plot(NRZ2);
title('Réponse en sortie du filtre émission sur la 1e chaine');

zr2_bis = filter(h2,1,NRZ2);

% tracé du diagramme de l'oeil
figure(4);
plot(reshape(zr2_bis,Ns,length(zr2_bis)/Ns));
title("Diagramme de l'oeil en sortie du filtre de la 1e chaine de transmission");


%% Introduction du bruit

TEB2 = zeros(1,9);
TEB2th =  zeros(1,9);

for i=0:8
    Px = mean(abs(NRZ2).^2);        % puissance du signal              
    SNRdb = i;
    SNR = 10^(SNRdb/10);       % SNR = EB/n0
    sigma = sqrt(Px*Ns/(2*log2(M)*SNR));

    bruit = sigma*randn(1,length(NRZ2));    % expression du bruit

    NRZ2_bruite = NRZ2 + bruit;
    
    zr2 = filter(h2,1,NRZ2_bruite);
    
    
    % signaux échantillonnés
    z2_ech = zr2(1,(n0_1:Ns:end));
    
    % décision
    d2 = sign(z2_ech);
    z_res_2 = (d2 + 1)/2;
    
    % calcul des taux d'erreur binaire
    nb_bits_errones_2 = length(find(z_res_2-bits ~= 0));
    nb_bits_totaux = length(bits);

    tx_erreur_binaire_2 = nb_bits_errones_2/nb_bits_totaux;
    
    TEB2(i+1) = tx_erreur_binaire_2;
    
    TEB2th(i+1) = qfunc(sqrt(SNR));
end


% tracé TEB théorique et expérimental en fonction du SNRdb

figure(5);
semilogy([0:8],TEB2,[0:8],TEB2th);
title('Evolution du TEB2 et du TEB2 théorique en fonction de SNRdb');

% on choisit initialement Nb = 10000 pour avoir une précision du TEB à 10%
% (voir annexe)



%% 4. Deuxième chaine

Ts3 = Tb*2;     % durée entre symboles
Ns3 = Fe*Ts3;  % nb d'échantillons utilisés par symbole
M = 4;  % ordre de la modulation
t0_3 = Ts3;
n0_3 = t0_3/Te;


h3 = ones(1,Ns3);

% Mapping
x3 = (2*bi2de(reshape(bits,2,length(bits)/2).')-3).';

% h et hr sont des rectangles (h3 = hr3)

somme3 = kron(x3,[1 zeros(1,Ns3-1)]);

NRZ3 = filter(h3,1,somme3);
figure(6);
plot(NRZ3);
title('Réponse en sortie du filtre émission sur la 2e chaine');

zr3_bis = filter(h3,1,NRZ3);

% tracé du diagramme de l'oeil
figure(7);
plot(reshape(zr3_bis,Ns3,length(zr3_bis)/Ns3));
title("Diagramme de l'oeil en sortie du filtre de la 2e chaine de transmission");



%% Introduction du bruit

TEB3 = zeros(1,9);
TEB3th =  zeros(1,9);

TES3 = zeros(1,9);
TES3th =  zeros(1,9);

for i=0:8
    Px = mean(abs(NRZ3).^2);        % puissance du signal              
    SNRdb = i;
    SNR = 10^(SNRdb/10);       % SNR = EB/n0
    sigma = sqrt(Px*Ns3/(2*log2(M)*SNR));

    bruit = sigma*randn(1,length(NRZ3));    % expression du bruit

    NRZ3_bruite = NRZ3 + bruit;
    
    zr3 = filter(h3,1,NRZ3_bruite);
    
    
    % signaux échantillonnés
    z3_ech = zr3(1,(n0_3:Ns3:end))/Ns3;
    
    % décision

    z3_ech(find(z3_ech <= -2)) = -3;
    z3_ech(find(z3_ech <= 0 & z3_ech > -2)) = -1;
    z3_ech(find(z3_ech <= 2 & z3_ech > 0)) = 1;
    z3_ech(find(z3_ech > 2)) = 3;

    SymbolesDecides = z3_ech;

    % Demapping
    z_res_3 = reshape(de2bi((SymbolesDecides + 3)/2).',1,length(bits));
    
    
    % calcul des taux d'erreur symbole
    nb_symboles = length(find(SymbolesDecides-x3 ~= 0));
    nb_symboles_totaux = length(x3);
    
    tx_erreur_symbole = nb_symboles/nb_symboles_totaux;
    
    TES3(i+1) = tx_erreur_symbole;
    TES3th(i+1) = 3/2*qfunc(sqrt((4/5)*SNR));
    

    % calcul des taux d'erreur binaire
    nb_bits_errones_3 = length(find(z_res_3-bits ~= 0));
    nb_bits_totaux = length(bits);

    tx_erreur_binaire_3 = nb_bits_errones_3/nb_bits_totaux;
    
    TEB3(i+1) = tx_erreur_binaire_3;
    TEB3th(i+1) = 3/4*qfunc(sqrt((4/5)*SNR));
end


% tracé TES théorique et expérimental en fonction du SNRdb

figure(8);
semilogy([0:8],TES3,[0:8],TES3th);
title('Evolution du TES3 et du TES3 théorique en fonction de SNRdb');


% tracé TEB théorique et expérimental en fonction du SNRdb

figure(9);
semilogy([0:8],TEB3,[0:8],TEB3th);
title('Evolution du TEB3 et du TEB3 théorique en fonction de SNRdb');

% on choisit initialement Nb = 10000 pour avoir une précision du TEB à 10%
% (voir annexe)


%% Tracé final
figure(10);
semilogy([0:8],TEB,'r',[0:8],TEBth,'-.r',[0:8],TEB2,'b',[0:8],TEB2th,'-.b',[0:8],TEB3,'g',[0:8],TEB3th,'-.g',[0:8],TES3,'k',[0:8],TES3th,'-.k');
title('Evolution des TEB/TES et des TEB/TES théorique en fonction de SNRdb');
legend('TEB1','TEB1th','TEB2','TEB2th','TEB3','TEB3th','TES3','TES3th');

% on choisit initialement Nb = 10000 pour avoir une précision du TEB à 10%
% (voir annexe)



%% Tracé pour efficacité spectrale

g1 = conv(h1,h1);
g2 = conv(h1,h2);
g3 = conv(h3,h3);

Hg1 = abs(fft(g1,N));
Hg2 = abs(fft(g2,N));
Hg3 = abs(fft(g3,N));

Frequences1 = 1/Te*(1:length(Hg1));
Frequences2 = 1/Te*(1:length(Hg2));
Frequences3 = 1/Te*(1:length(Hg3));

figure(11);
plot(Frequences1,Hg1/max(Hg1),Frequences2,Hg2/max(Hg2),Frequences3,Hg3/max(Hg3));
title('Comparaison des densités spectrales des différentes chaines');
legend('Hg1/max(Hg1)','Hg2/max(Hg2)','Hg3/max(Hg3)');
xlabel('Fréquences(Hz)');
ylabel('fft');
