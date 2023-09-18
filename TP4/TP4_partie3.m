clear
close all

Fe = 10000;  % fréquence d'échantillonnage
Rb = 2000;   % débit binaire 
Tb = 1/Rb;
Nb = 10000;
Te = 1/Fe;

bits = randi([0,1],1,Nb);  % Génération de Nb bits


M = 4;  % ordre de la modulation car chaine QPSK => 4 symboles
Ts = log2(M)/Rb;     % durée entre symboles
Ns = Fe*Ts;  % nb d'échantillons utilisés par symbole
fp = 2000; % fréquence porteuse


t0_1 = Ts;
n0_1 = 1;  %t0_1/Te;  => on enlève le décalage


%% 3.1 Implantation de la chaine sur fréquence porteuse (chaine de transmission QPSK)

alpha = 0.35;

% Mapping de Gray
x1 = bi2de(reshape(bits,Nb/2,2))';
ak = zeros(1,length(x1));
bk = zeros(1,length(x1));

% Génération des symboles ak
ak(find(x1 == 0 | x1 == 1)) = -1;
ak(find(x1 == 2 | x1 == 3)) = 1;

% Génération des symboles bk
bk(find(x1 == 0 | x1 == 2)) = -1;
bk(find(x1 == 1 | x1 == 3)) = 1;

dk = ak + 1i*bk;


% réponse impulsionnelle
h1 = rcosdesign(alpha,8,Ns);
N = length(h1);

% Suréchantillonnage
somme1 = kron(dk,[1 zeros(1,Ns-1)]);

% Enveloppe complexe
xe_1 = filter(h1,1,[somme1 zeros(1,(N-1)/2)]);
xe_1 = xe_1((N-1)/2+1:end);

I = real(xe_1);
Q = imag(xe_1);

tps = [0:Te:(length(xe_1)-1)*Te];
x_1 = real(xe_1.*exp(2*1i*pi*fp*tps));


%% 3.1.1) Tracé des signaux
figure();
plot(tps,I);
title('Signal généré sur la voie en phase');
xlabel('Temps (en s)');
ylabel('I(t)');

figure();
plot(tps,Q);
title('Signal généré sur la voie en quadrature');
xlabel('Temps (en s)');
ylabel('Q(t)');

figure();
plot(tps,x_1);
title('Signal généré sur fréquence porteuse');
xlabel('Temps (en s)');
ylabel('x_1(t)');



%% 3.1.2) Estimer et tracer la densité spectrale de puissance du signal modulé sur fréquence porteuse
dsp1 = pwelch(x_1,[],[],[],Fe,'twosided');
Frequences = (Fe/length(dsp1))*[0:1:length(dsp1)-1];

figure();
plot(Frequences, 10*log10(dsp1));
title('DSP du signal modulé sur fréquence porteuse');
xlabel('Fréquences (en Hz)');
ylabel('Puissance/Frequence (en dB/Hz)');
grid on;



%% 3.1.3) Implantation chaîne complète sans bruit

% Retour en bande de base

cosx_1 = 2*cos(2*pi*fp*tps).*x_1;
sinx_1 = 2*sin(2*pi*fp*tps).*x_1;
x1_ret = cosx_1 - 1i*sinx_1;


% Démodulation bande de base
zr1 = filter(h1,1,[x1_ret zeros(1,(N-1)/2)]);
zr1 = zr1((N-1)/2+1:end);

% tracé du diagramme de l'oeil
figure();
plot(reshape(real(zr1),Ns,length(real(zr1))/Ns));
title("Diagramme de l'oeil en sortie du filtre de la chaine de transmission");

% % signal échantillonné
% z1_ech = zr1(1,(n0_1:Ns:end));
%     
% % décision
% d1 = angle(z1_ech);
% decision(find(d1 >= 0 & d1 <= pi/2)) = 3;
% decision(find(d1 >= -pi & d1 <= -pi/2)) = 0;
% decision(find(d1 >= pi/2 & d1 <= pi)) = 1;
% decision(find(d1 >= -pi/2 & d1 <= 0)) = 2;
% 
% 
% % demapping
% z_res_1 = reshape(de2bi(decision),1,length(bits));
% 
% % calcul des taux d'erreur binaire
% nb_bits_errones_1 = length(find(z_res_1-bits ~= 0));
% nb_bits_totaux = length(bits);
% 
% tx_erreur_binaire_1 = nb_bits_errones_1/nb_bits_totaux




%% 3.1.4) Introduction du bruit blanc et TEB

TEB = zeros(1,9);
TEBth =  zeros(1,9);


for i=0:8
    
    bits = randi([0,1],1,Nb);   % on regénère les bits à chaque itération

    % Mapping de Gray
    x1 = bi2de(reshape(bits,Nb/2,2))';
    ak = zeros(1,length(x1));
    bk = zeros(1,length(x1));

    % Génération des symboles ak
    ak(find(x1 == 0 | x1 == 1)) = -1;
    ak(find(x1 == 2 | x1 == 3)) = 1;

    % Génération des symboles bk
    bk(find(x1 == 0 | x1 == 2)) = -1;
    bk(find(x1 == 1 | x1 == 3)) = 1;

    dk = ak + 1i*bk;


    % réponse impulsionnelle
    h1 = rcosdesign(alpha,8,Ns);
    N = length(h1);

    % Suréchantillonnage
    somme1 = kron(dk,[1 zeros(1,Ns-1)]);

    % Enveloppe complexe
    xe_1 = filter(h1,1,[somme1 zeros(1,(N-1)/2)]);
    xe_1 = xe_1((N-1)/2+1:end);

    I = real(xe_1);

    tps = [0:Te:(length(xe_1)-1)*Te];
    x_1 = real(xe_1.*exp(2*1i*pi*fp*tps));

    
    Px = mean(abs(x_1).^2);        % puissance du signal              
    SNRdb = i;
    SNR = 10^(SNRdb/10);       % SNR = EB/n0
    sigma = sqrt(Px*Ns/(2*log2(M)*SNR));

    bruit = sigma*randn(1,length(x_1));    % expression du bruit

    xe1_bruite = x_1 + bruit;
    
    
    % Retour en bande de base

    cosx_1 = 2*cos(2*pi*fp*tps).*xe1_bruite;
    sinx_1 = 2*sin(2*pi*fp*tps).*xe1_bruite;
    x1_ret = cosx_1 - 1i*sinx_1;


    % Démodulation bande de base
    zr1 = filter(h1,1,[x1_ret zeros(1,(N-1)/2)]);
    zr1 = zr1((N-1)/2+1:end);
    
    
    % signaux échantillonnés
    z1_ech = zr1(1,(n0_1:Ns:end));

    % décision
    d1 = angle(z1_ech);
    decision(find(d1 >= 0 & d1 <= pi/2)) = 3;
    decision(find(d1 >= -pi & d1 <= -pi/2)) = 0;
    decision(find(d1 >= pi/2 & d1 <= pi)) = 1;
    decision(d1 >= -pi/2 & d1 <= 0) = 2;

    % demapping
    z_res_1 = reshape(de2bi(decision),1,length(bits));

    % calcul des taux d'erreur binaire
    nb_bits_errones_1 = length(find(z_res_1-bits ~= 0));
    nb_bits_totaux = length(bits);

    tx_erreur_binaire_1 = nb_bits_errones_1/nb_bits_totaux;
    
    TEB(i+1) = tx_erreur_binaire_1;
    
    TEBth(i+1) = qfunc(sqrt(2*SNR*log2(M))*sin(pi/M));
end


%% 3.1.5) Tracé TEB théorique et expérimental en fonction du SNRdb

figure();
semilogy([0:8],TEB,[0:8],TEBth);
title('Evolution du TEB et du TEB théorique en fonction de SNRdb');






%% 3.2 Implantation de la chaine passe-bas ́equivalente

alpha = 0.35;

% Mapping de Gray
x2 = bi2de(reshape(bits,Nb/2,2))';
ak = zeros(1,length(x1));
bk = zeros(1,length(x1));

% Génération des symboles ak
ak(find(x2 == 0 | x2 == 1)) = -1;
ak(find(x2 == 2 | x2 == 3)) = 1;

% Génération des symboles bk
bk(find(x2 == 0 | x2 == 2)) = -1;
bk(find(x2 == 1 | x2 == 3)) = 1;

dk = ak + 1i*bk;

% même réponse impulsionnelle et suréchantillonnage

% Enveloppe complexe
x_2 = filter(h1,1,[somme1 zeros(1,(N-1)/2)]);
x_2 = x_2((N-1)/2+1:end);

I = real(x_2);
Q = imag(x_2);

tps = [0:Te:(length(x_2)-1)*Te];


%% 3.2.1) Tracé des signaux

figure();
plot(tps,I);
title('Signal généré sur la voie en phase');
xlabel('Temps (en s)');
ylabel('I(t)');

figure();
plot(tps,Q);
title('Signal généré sur la voie en quadrature');
xlabel('Temps (en s)');
ylabel('Q(t)');

figure();
plot(tps,real(x_2));
title('Signal généré sur fréquence porteuse');
xlabel('Temps (en s)');
ylabel('x_2(t)');



%% 3.2.2) Estimer et tracer la densité spectrale de puissance du signal modulé sur fréquence porteuse

dsp2 = pwelch(x_2,[],[],[],Fe,'twosided');
Frequences = (Fe/length(dsp2))*[0:1:length(dsp2)-1];

figure();
plot(Frequences, 10*log(dsp2));
title('DSP du signal modulé sur fréquence porteuse');
xlabel('Fréquences (en Hz)');
ylabel('Puissance/Frequence (en dB/Hz)');
grid on;


%% 3.2.3) Implantation chaîne complète sans bruit

% Démodulation bande de base
zr2 = filter(h1,1,[x_2 zeros(1,(N-1)/2)]);
zr2 = zr2((N-1)/2+1:end);

% tracé du diagramme de l'oeil
figure();
plot(reshape(real(zr2),Ns,length(real(zr2))/Ns));
title("Diagramme de l'oeil en sortie du filtre de la chaine de transmission");

% 
% % signal échantillonné
% z2_ech = zr2(1,(n0_1:Ns:end));
%     
% % décision
% d2 = angle(z2_ech);
% decision2(find(d2 >= 0 & d2 <= pi/2)) = 3;
% decision2(find(d2 >= -pi & d2 <= -pi/2)) = 0;
% decision2(find(d2 >= pi/2 & d2 <= pi)) = 1;
% decision2(find(d2 >= -pi/2 & d2 <= 0)) = 2;
% 
% 
% % demapping
% z_res_2 = reshape(de2bi(decision2),1,length(bits));
% 
% % calcul des taux d'erreur binaire
% nb_bits_errones_2 = length(find(z_res_2-bits ~= 0));
% nb_bits_totaux = length(bits);
% 
% tx_erreur_binaire_2 = nb_bits_errones_2/nb_bits_totaux




%% 3.2.4) Introduction du bruit blanc et TEB

TEB2 = zeros(1,9);
TEBth =  zeros(1,9);


for i=0:8
    
    bits = randi([0,1],1,Nb);   % on regénère les bits à chaque itération
    
    % Mapping de Gray
    x2 = bi2de(reshape(bits,Nb/2,2))';
    ak = zeros(1,length(x2));
    bk = zeros(1,length(x2));

    % Génération des symboles ak
    ak(find(x2 == 0 | x2 == 1)) = -1;
    ak(find(x2 == 2 | x2 == 3)) = 1;

    % Génération des symboles bk
    bk(find(x2 == 0 | x2 == 2)) = -1;
    bk(find(x2 == 1 | x2 == 3)) = 1;

    dk = ak + 1i*bk;
    
    if (i==1)
        figure();
        plot(ak,bk,'*');
        title("Tracé de la constellation en sortie du mapping");
    end
    
    
    % même réponse impulsionnelle
    
    % Suréchantillonnage
    somme2 = kron(dk,[1 zeros(1,Ns-1)]);

    % Enveloppe complexe
    x_2 = filter(h1,1,[somme2 zeros(1,(N-1)/2)]);
    x_2 = x_2((N-1)/2+1:end);

    I = real(x_2);
    Q = imag(x_2);

    tps = [0:Te:(length(x_2)-1)*Te];

    
    Px = mean(abs(x_2).^2);        % puissance du signal              
    SNRdb = i;
    SNR = 10^(SNRdb/10);       % SNR = EB/n0
    sigma = sqrt(Px*Ns/(2*log2(M)*SNR));

    bruit_I = sigma*randn(1,length(x_2));    % expression du bruit
    bruit_Q = sigma*randn(1,length(x_2));
    bruit = bruit_I +1i*bruit_Q;
    
    xe2_bruite = x_2 + bruit;

    % Démodulation bande de base
    zr2 = filter(h1,1,[xe2_bruite zeros(1,(N-1)/2)]);
    zr2 = zr2((N-1)/2+1:end);
    
    % signaux échantillonnés
    z2_ech = zr2(1,(n0_1:Ns:end));

    
    % décision
    d2 = angle(z2_ech);
    decision2(find(d2 >= 0 & d2 <= pi/2)) = 3;
    decision2(find(d2 >= -pi & d2 <= -pi/2)) = 0;
    decision2(find(d2 >= pi/2 & d2 <= pi)) = 1;
    decision2(find(d2 >= -pi/2 & d2 <= 0)) = 2;

    if (i==1)
        figure();
        plot(real(z2_ech),imag(z2_ech),'*');
        title("Tracé de la constellation en sortie de l'échantillonneur");
    end
    
    % demapping
    z_res_2 = reshape(de2bi(decision2),1,length(bits));

    % calcul des taux d'erreur binaire
    nb_bits_errones_2 = length(find(z_res_2-bits ~= 0));
    nb_bits_totaux = length(bits);

    tx_erreur_binaire_2 = nb_bits_errones_2/nb_bits_totaux;
    
    TEB2(i+1) = tx_erreur_binaire_2;
    
    TEBth(i+1) = qfunc(sqrt(2*SNR*log2(M))*sin(pi/M));
end


% Tracé TEB théorique et expérimental en fonction du SNRdb

figure();
semilogy([0:8],TEB2,[0:8],TEBth);
title('Evolution du TEB et du TEB théorique en fonction de SNRdb');


%% 3.2.5 Tracé des constellations

% (voir figure plus haut)

%% 3.2.6 Tracé TEB sur fréquence porteuse et chaine passe bas équivalente

figure();
plot([0:8],TEB,'b',[0:8],TEB2,'r');
title('Evolution des 2 TEB expérimentaux calculés en fonction de SNRdb');
xlabel('SNRdb');
ylabel('TEB');
legend('Chaine sur fréquence porteuse','Chaine passe bas équivalente');
grid on;



