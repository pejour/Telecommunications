clear
close all

Fe = 12000;  % fréquence d'échantillonnage
Rb = 6000;   % débit binaire 
Tb = 1/Rb;
Nb = 12000;
Te = 1/Fe;
bits = randi([0,1],1,Nb);  % Génération de Nb bits
fp = 2000; % fréquence porteuse




%% 4.1 Etude de chaque chaine de transmission


%% Modulation 4-ASK

alpha = 0.35;
M1 = 4;       % ordre de la modulation car chaine 4-ASK => 4 symboles
Ts1 = log2(M1)/Rb;    % durée entre symboles
Ns1 = Ts1*Fe;    % nb d'échantillons utilisés par symbole

t0_1 = Ts1;
n0_1 = 1;  %t0_1/Te;  => on enlève le décalage

% Filtres d'émission et de réception
h1 = rcosdesign(alpha,8,Ns1);
N = length(h1);

% Mapping
symboles_4_ASK = mapping_4_ASK(bits);


% Suréchantillonnage
somme_4_ASK = kron(symboles_4_ASK,[1 zeros(1,Ns1-1)]);

% Enveloppe complexe
xe_1 = filter(h1,1,[somme_4_ASK zeros(1,(N-1)/2)]);
xe_1 = xe_1((N-1)/2+1:end);

tps = [0:Te:(length(xe_1)-1)*Te];
x_1 = real(xe_1.*exp(2*1i*pi*fp*tps));


% 1) Implantation chaîne complète sans bruit

cosx_1 = 2*cos(2*pi*fp*tps).*x_1;
sinx_1 = 2*sin(2*pi*fp*tps).*x_1;
x1_ret = cosx_1 - 1i*sinx_1;


% Démodulation bande de base
zr1 = filter(h1,1,[x1_ret zeros(1,(N-1)/2)]);
zr1 = zr1((N-1)/2+1:end);


% tracé du diagramme de l'oeil
figure();
plot(reshape(real(zr1),Ns1,length(real(zr1))/Ns1));
title("Diagramme de l'oeil en sortie du filtre de la chaine de transmission (4-ASK)");


% % signal échantillonné
% z1_ech = zr1(1,(n0_1:Ns1:end));
%     
% % décision
% decision1 = decision_4_ASK(z1_ech);
% 
% % demapping
% z_res_1 = demapping_4_ASK(decision1);
% 
% % calcul des taux d'erreur binaire
% nb_bits_errones_1 = length(find(z_res_1-bits ~= 0));
% nb_bits_totaux = length(bits);
% 
% tx_erreur_binaire_1 = nb_bits_errones_1/nb_bits_totaux
 

 
% 2) Introduction du bruit blanc et TEB

TEB1 = zeros(1,9);
TEB1th =  zeros(1,9);


for i=0:8
    
    bits = randi([0,1],1,Nb);   % on regénère les bits à chaque itération
    
    % Mapping
    symboles_4_ASK = mapping_4_ASK(bits);
    
    % Tracé de la constellation en sortie du mapping
    if (i==1)
        figure();
        plot(symboles_4_ASK, zeros(1, length(symboles_4_ASK)), ['*' 0.5*(rand(1,3)+1)]);
        title('Constellation en sortie du mapping pour une modulation 4-ASK')
        xlabel('Symboles ak');
        ylabel('Symboles bk');
        grid on;
    end
    
    % Suréchantillonnage
    somme_4_ASK = kron(symboles_4_ASK,[1 zeros(1,Ns1-1)]);

    % Enveloppe complexe
    xe_1 = filter(h1,1,[somme_4_ASK zeros(1,(N-1)/2)]);
    xe_1 = xe_1((N-1)/2+1:end);

    tps = [0:Te:(length(xe_1)-1)*Te];
    x_1 = real(xe_1.*exp(2*1i*pi*fp*tps));

    
    Px = mean(abs(x_1).^2);        % puissance du signal              
    SNRdb = i;
    SNR = 10^(SNRdb/10);       % SNR = EB/n0
    sigma = sqrt(Px*Ns1/(2*log2(M1)*SNR));

    bruit = sigma*randn(1,length(x_1));    % expression du bruit
    
    xe1_bruite = x_1 + bruit;

    % Retour en bande de base

    cosx_1 = 2*cos(2*pi*fp*tps).*xe1_bruite;
    sinx_1 = 2*sin(2*pi*fp*tps).*xe1_bruite;
    x1_ret = cosx_1 - 1i*sinx_1;


    % Démodulation bande de base
    zr1 = filter(h1,1,[x1_ret zeros(1,(N-1)/2)]);
    zr1 = zr1((N-1)/2+1:end);
    
    % signal échantillonné
    z1_ech = zr1(1,(n0_1:Ns1:end));

    % décision
    decision1 = decision_4_ASK(z1_ech);

    % Tracé de la constellation en sortie de l'échantillonneur
    if (i==1)
        figure();
        plot(real(z1_ech), zeros(1, length(z1_ech)), ['*' 0.5*(rand(1,3)+1)]);
        title("Constellation en sortie de l'échantillonneur pour une modulation 4-ASK")
        xlabel('Symboles ak');
        ylabel('Symboles bk');
        grid on;
    end
    
    % demapping
    z_res_1 = demapping_4_ASK(decision1);
    

    % calcul des taux d'erreur binaire
    nb_bits_errones_1 = length(find(z_res_1-bits ~= 0));
    nb_bits_totaux = length(bits);

    tx_erreur_binaire_1 = nb_bits_errones_1/nb_bits_totaux;
    
    TEB1(i+1) = tx_erreur_binaire_1;
    
    TEB1th(i+1) = ((M1-1)/M1) * qfunc(sqrt((6*log2(M1)*SNR)/(M1*M1-1)));
end


% Tracé TEB théorique et expérimental en fonction du SNRdb

figure();
semilogy([0:8],TEB1,[0:8],TEB1th);
title('Evolution du TEB1 et du TEB1 théorique en fonction de SNRdb');




%% Modulation QPSK

alpha = 0.35;
M2 = 4;       % ordre de la modulation car chaine QPSK => 4 symboles
Ts2 = log2(M2)/Rb;    % durée entre symboles
Ns2 = Ts2*Fe;    % nb d'échantillons utilisés par symbole

t0_2 = Ts2;
n0_2 = 1;  %t0_2/Te;  => on enlève le décalage

% Filtres d'émission et de réception
h2 = rcosdesign(alpha,8,Ns2);
N = length(h2);

% Mapping + Tracé de la constellation en sortie du mapping
symboles_QPSK = transpose(qammod(transpose(bits), M2,'gray', 'InputType', 'bit','PlotConstellation',true));

% Suréchantillonnage
somme_QPSK = kron(symboles_QPSK,[1 zeros(1,Ns2-1)]);

% Enveloppe complexe
xe_2 = filter(h2,1,[somme_QPSK zeros(1,(N-1)/2)]);
xe_2 = xe_2((N-1)/2+1:end);

tps = [0:Te:(length(xe_2)-1)*Te];
x_2 = real(xe_2.*exp(2*1i*pi*fp*tps));


% 1) Implantation chaîne complète sans bruit

cosx_2 = 2*cos(2*pi*fp*tps).*x_2;
sinx_2 = 2*sin(2*pi*fp*tps).*x_2;
x2_ret = cosx_2 - 1i*sinx_2;


% Démodulation bande de base
zr2 = filter(h2,1,[x2_ret zeros(1,(N-1)/2)]);
zr2 = zr2((N-1)/2+1:end);


% tracé du diagramme de l'oeil
figure();
plot(reshape(real(zr2),Ns2,length(real(zr2))/Ns2));
title("Diagramme de l'oeil en sortie du filtre de la chaine de transmission (QPSK)");


% % signal échantillonné
% z2_ech = zr2(1,(n0_2:Ns2:end));
% 
% % demapping
% z_res_2 = transpose(qamdemod(transpose(z2_ech), M2));
% 
% % calcul des taux d'erreur binaire
% nb_bits_errones_2 = length(find(z_res_2-bits ~= 0));
% nb_bits_totaux = length(bits);
% 
% tx_erreur_binaire_2 = nb_bits_errones_2/nb_bits_totaux;
%  

 
% 2) Introduction du bruit blanc et TEB

TEB2 = zeros(1,9);
TEB2th =  zeros(1,9);


for i=0:8
    
    bits = randi([0,1],1,Nb);   % on regénère les bits à chaque itération
    
    % Mapping
    symboles_QPSK = transpose(qammod(transpose(bits), M2));
    
    
    % Suréchantillonnage
    somme_QPSK = kron(symboles_QPSK,[1 zeros(1,Ns2-1)]);

    % Enveloppe complexe
    xe_2 = filter(h2,1,[somme_QPSK zeros(1,(N-1)/2)]);
    xe_2 = xe_2((N-1)/2+1:end);

    tps = [0:Te:(length(xe_2)-1)*Te];
    x_2 = real(xe_2.*exp(2*1i*pi*fp*tps));

    
    Px = mean(abs(x_2).^2);        % puissance du signal              
    SNRdb = i;
    SNR = 10^(SNRdb/10);       % SNR = EB/n0
    sigma = sqrt(Px*Ns2/(2*log2(M2)*SNR));

    bruit = sigma*randn(1,length(x_2));    % expression du bruit
    
    xe2_bruite = x_2 + bruit;

    % Retour en bande de base

    cosx_2 = 2*cos(2*pi*fp*tps).*xe2_bruite;
    sinx_2 = 2*sin(2*pi*fp*tps).*xe2_bruite;
    x2_ret = cosx_2 - 1i*sinx_2;


    % Démodulation bande de base
    zr2 = filter(h2,1,[x2_ret zeros(1,(N-1)/2)]);
    zr2 = zr2((N-1)/2+1:end);
    
    % signal échantillonné
    z2_ech = zr2(1,(n0_2:Ns2:end));

    % Tracé de la constellation en sortie de l'échantillonneur
    if (i==1)
        figure();
        plot(z2_ech, ['*' 0.5*(rand(1,3)+1)],'MarkerSize',5,'LineWidth',2);
        title("Constellation en sortie de l'échantillonneur pour une modulation QPSK")
        xlabel('Symboles ak');
        ylabel('Symboles bk');
        grid on;
    end
    
    % demapping
    z_res_2 = transpose(qamdemod(transpose(z2_ech), M2));

    % calcul des taux d'erreur binaire
    nb_bits_errones_2 = length(find(z_res_2-bits ~= 0));
    nb_bits_totaux = length(bits);

    tx_erreur_binaire_2 = nb_bits_errones_2/nb_bits_totaux;
    
    TEB2(i+1) = tx_erreur_binaire_2;
    
    TEB2th(i+1) = 2*(1-1/sqrt(M2)) * qfunc(sqrt((3*log2(M2)*SNR)/(M2-1)));
end


% Tracé TEB théorique et expérimental en fonction du SNRdb

figure();
semilogy([0:8],TEB2,[0:8],TEB2th);
title('Evolution du TEB2 et du TEB2 théorique en fonction de SNRdb');








%% Modulation 8-PSK

alpha = 0.35;
M3 = 8;       % ordre de la modulation car chaine 8-ASK => 8 symboles
Ts3 = log2(M3)/Rb;    % durée entre symboles
Ns3 = Ts3*Fe;    % nb d'échantillons utilisés par symbole

t0_3 = Ts3;
n0_3 = 1;  %t0_3/Te;  => on enlève le décalage

% Filtres d'émission et de réception
h3 = rcosdesign(alpha,8,Ns3);
N = length(h3);

% Mapping
symboles_8_PSK = mapping_8_PSK(bits);


% Suréchantillonnage
somme_8_PSK = kron(symboles_8_PSK,[1 zeros(1,Ns3-1)]);


% Enveloppe complexe
xe_3 = filter(h3,1,[somme_8_PSK zeros(1,(N-1)/2)]);
xe_3 = xe_3((N-1)/2+1:end);

tps = [0:Te:(length(xe_3)-1)*Te];
x_3 = real(xe_3.*exp(2*1i*pi*fp*tps));


% 1) Implantation chaîne complète sans bruit

cosx_3 = 2*cos(2*pi*fp*tps).*x_3;
sinx_3 = 2*sin(2*pi*fp*tps).*x_3;
x3_ret = cosx_3 - 1i*sinx_3;


% Démodulation bande de base
zr3 = filter(h3,1,[x3_ret zeros(1,(N-1)/2)]);
zr3 = zr3((N-1)/2+1:end);


% tracé du diagramme de l'oeil
figure();
plot(reshape(real(zr3),Ns3,length(real(zr3))/Ns3));
title("Diagramme de l'oeil en sortie du filtre de la chaine de transmission (8-PSK)");


% % signal échantillonné
% z3_ech = zr3(1,(n0_3:Ns3:end));
%     
% % décision
% decision3 = decision_8_PSK(z3_ech);
% 
% % demapping
% z_res_3 = demapping_8_PSK(decision3);
% 
% % calcul des taux d'erreur binaire
% nb_bits_errones_3 = length(find(z_res_3-bits ~= 0));
% nb_bits_totaux = length(bits);
% 
% tx_erreur_binaire_3 = nb_bits_errones_3/nb_bits_totaux
 

 
% 2) Introduction du bruit blanc et TEB

TEB3 = zeros(1,9);
TEB3th =  zeros(1,9);


for i=0:8
    
    bits = randi([0,1],1,Nb);   % on regénère les bits à chaque itération
    
    % Mapping
    symboles_8_PSK = mapping_8_PSK(bits);
    
    % Tracé de la constellation en sortie du mapping
    if (i==1)
        figure();
        plot(symboles_8_PSK, ['*' 0.5*(rand(1,3)+1)],'MarkerSize',10,'LineWidth',2);
        title('Constellation en sortie du mapping pour une modulation 8-PSK')
        xlabel('Symboles ak');
        ylabel('Symboles bk');
        grid on;
    end
    
    % Suréchantillonnage
    somme_8_PSK = kron(symboles_8_PSK,[1 zeros(1,Ns3-1)]);

    % Enveloppe complexe
    xe_3 = filter(h3,1,[somme_8_PSK zeros(1,(N-1)/2)]);
    xe_3 = xe_3((N-1)/2+1:end);

    tps = [0:Te:(length(xe_3)-1)*Te];
    x_3 = real(xe_3.*exp(2*1i*pi*fp*tps));

    
    Px = mean(abs(x_3).^2);        % puissance du signal              
    SNRdb = i;
    SNR = 10^(SNRdb/10);       % SNR = EB/n0
    sigma = sqrt(Px*Ns3/(2*log2(M3)*SNR));

    bruit = sigma*randn(1,length(x_3));    % expression du bruit
    
    xe3_bruite = x_3 + bruit;

    % Retour en bande de base

    cosx_3 = 2*cos(2*pi*fp*tps).*xe3_bruite;
    sinx_3 = 2*sin(2*pi*fp*tps).*xe3_bruite;
    x3_ret = cosx_3 - 1i*sinx_3;


    % Démodulation bande de base
    zr3 = filter(h3,1,[x3_ret zeros(1,(N-1)/2)]);
    zr3 = zr3((N-1)/2+1:end);
    
    % signal échantillonné
    z3_ech = zr3(1,(n0_1:Ns3:end));

    % décision
    decision3 = decision_8_PSK(z3_ech);

    % Tracé de la constellation en sortie de l'échantillonneur
    if (i==1)
        figure();
        plot(z3_ech, ['*' 0.5*(rand(1,3)+1)],'MarkerSize',5,'LineWidth',2);
        title("Constellation en sortie de l'échantillonneur pour une modulation 8-PSK")
        xlabel('Symboles ak');
        ylabel('Symboles bk');
        grid on;
    end
    
    % demapping
    z_res_3 = demapping_8_PSK(decision3);
    

    % calcul des taux d'erreur binaire
    nb_bits_errones_3 = length(find(z_res_3-bits ~= 0));
    nb_bits_totaux = length(bits);

    tx_erreur_binaire_3 = nb_bits_errones_3/nb_bits_totaux;
    
    TEB3(i+1) = tx_erreur_binaire_3;
    
    TEB3th(i+1) = (2/3)*qfunc(sqrt(6*SNR)*sin(pi/M3));
end


% Tracé TEB théorique et expérimental en fonction du SNRdb

figure();
semilogy([0:8],TEB3,[0:8],TEB3th);
title('Evolution du TEB3 et du TEB3 théorique en fonction de SNRdb');





%% Modulation 16-QAM

alpha = 0.35;
M4 = 16;       % ordre de la modulation car chaine 16-QAM => 16 symboles
Ts4 = log2(M4)/Rb;    % durée entre symboles
Ns4 = Ts4*Fe;    % nb d'échantillons utilisés par symbole

t0_4 = Ts4;
n0_4 = 1;  %t0_4/Te;  => on enlève le décalage

% Filtres d'émission et de réception
h4 = rcosdesign(alpha,8,Ns4);
N = length(h4);

% Mapping + Tracé de la constellation en sortie du mapping
symboles_16_QAM = transpose(qammod(transpose(bits), M4,'gray', 'InputType', 'bit','PlotConstellation',true));

% Suréchantillonnage
somme_16_QAM = kron(symboles_16_QAM,[1 zeros(1,Ns4-1)]);

% Enveloppe complexe
xe_4 = filter(h4,1,[somme_16_QAM zeros(1,(N-1)/2)]);
xe_4 = xe_4((N-1)/2+1:end);

tps = [0:Te:(length(xe_4)-1)*Te];
x_4 = real(xe_4.*exp(2*1i*pi*fp*tps));


% 1) Implantation chaîne complète sans bruit

cosx_4 = 2*cos(2*pi*fp*tps).*x_4;
sinx_4 = 2*sin(2*pi*fp*tps).*x_4;
x4_ret = cosx_4 - 1i*sinx_4;


% Démodulation bande de base
zr4 = filter(h4,1,[x4_ret zeros(1,(N-1)/2)]);
zr4 = zr4((N-1)/2+1:end);


% tracé du diagramme de l'oeil
figure();
plot(reshape(real(zr4),Ns4,length(real(zr4))/Ns4));
title("Diagramme de l'oeil en sortie du filtre de la chaine de transmission (16-QAM)");


% % signal échantillonné
% z4_ech = zr4(1,(n0_4:Ns4:end));
% 
% % demapping
% z_res_4 = transpose(qamdemod(transpose(z4_ech), M4, 'gray', 'OutputType', 'bit'));
% 
% % calcul des taux d'erreur binaire
% nb_bits_errones_4 = length(find(z_res_4-bits ~= 0));
% nb_bits_totaux = length(bits);
% 
% tx_erreur_binaire_4 = nb_bits_errones_4/nb_bits_totaux


 
% 2) Introduction du bruit blanc et TEB

TEB4 = zeros(1,9);
TEB4th =  zeros(1,9);


for i=0:8
    
    bits = randi([0,1],1,Nb);   % on regénère les bits à chaque itération
    
    % Mapping
    symboles_16_QAM = transpose(qammod(transpose(bits), M4));
    
    
    % Suréchantillonnage
    somme_16_QAM = kron(symboles_16_QAM,[1 zeros(1,Ns4-1)]);

    % Enveloppe complexe
    xe_4 = filter(h4,1,[somme_16_QAM zeros(1,(N-1)/2)]);
    xe_4 = xe_4((N-1)/2+1:end);

    tps = [0:Te:(length(xe_4)-1)*Te];
    x_4 = real(xe_4.*exp(2*1i*pi*fp*tps));

    
    Px = mean(abs(x_4).^2);        % puissance du signal              
    SNRdb = i;
    SNR = 10^(SNRdb/10);       % SNR = EB/n0
    sigma = sqrt(Px*Ns4/(2*log2(M4)*SNR));

    bruit = sigma*randn(1,length(x_4));    % expression du bruit
    
    xe4_bruite = x_4 + bruit;

    % Retour en bande de base

    cosx_4 = 2*cos(2*pi*fp*tps).*xe4_bruite;
    sinx_4 = 2*sin(2*pi*fp*tps).*xe4_bruite;
    x4_ret = cosx_4 - 1i*sinx_4;


    % Démodulation bande de base
    zr4 = filter(h4,1,[x4_ret zeros(1,(N-1)/2)]);
    zr4 = zr4((N-1)/2+1:end);
    
    % signal échantillonné
    z4_ech = zr4(1,(n0_4:Ns4:end));

    % Tracé de la constellation en sortie de l'échantillonneur
    if (i==1)
        figure();
        plot(z4_ech, ['*' 0.5*(rand(1,3)+1)],'MarkerSize',5,'LineWidth',2);
        title("Constellation en sortie de l'échantillonneur pour une modulation 16-QAM")
        xlabel('Symboles ak');
        ylabel('Symboles bk');
        grid on;
    end
    
    % demapping
    z_res_4 = transpose(qamdemod(transpose(z4_ech), M4));

    % calcul des taux d'erreur binaire
    nb_bits_errones_4 = length(find(z_res_4-bits ~= 0));
    nb_bits_totaux = length(bits);

    tx_erreur_binaire_4 = nb_bits_errones_4/nb_bits_totaux;
    
    TEB4(i+1) = tx_erreur_binaire_4;
    
    TEB4th(i+1) = (1-1/sqrt(M4))*qfunc(sqrt((3*log2(M4)*SNR)/(M4-1)));
end


% Tracé TEB théorique et expérimental en fonction du SNRdb

figure();
semilogy([0:8],TEB4,[0:8],TEB4th);
title('Evolution du TEB4 et du TEB4 théorique en fonction de SNRdb');






%% Comparaison des différents TEB calculés pour chaque modulation

figure();
plot([0:8],TEB1th,'r',[0:8],TEB2th,'b',[0:8],TEB3th,'g',[0:8],TEB4th,'m');
title('Evolution des différents TEB calculés en fonction de SNRdb');
xlabel('SNR en dB');
ylabel('TEB');
legend('4-ASK','QPSK','8-PSK','16-QAM');



%% Comparaison des DSP des signaux émis

dsp1 = pwelch(real(x_1),[],[],[],Fe,'twosided');
dsp2 = pwelch(real(x_2),[],[],[],Fe,'twosided');
dsp3 = pwelch(real(x_3),[],[],[],Fe,'twosided');
dsp4 = pwelch(real(x_4),[],[],[],Fe,'twosided');

Frequences1 = (Fe/length(dsp1))*[0:1:length(dsp1)-1];
Frequences2 = (Fe/length(dsp2))*[0:1:length(dsp2)-1];
Frequences3 = (Fe/length(dsp3))*[0:1:length(dsp3)-1];
Frequences4 = (Fe/length(dsp4))*[0:1:length(dsp4)-1];

figure();
plot(Frequences1, 10*log10(dsp1),'r', Frequences2, 10*log10(dsp2),'b', Frequences3, 10*log10(dsp3),'g', Frequences4, 10*log10(dsp4),'m');
title('Comparaison des DSP des signaux émis');
xlabel('Fréquences (en Hz)');
ylabel('Puissance/Frequence (en dB/Hz)');
legend('4-ASK','QPSK','8-PSK','16-QAM');
grid on;


