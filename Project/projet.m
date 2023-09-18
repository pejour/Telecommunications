close all;
clear;

% Création des figures
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);
d = 30;
figure('Name',"Signaux générés",'Position',[0,H/4,2*L/5,H/2]);
figure('Name',"DSP des signaux",'Position',[d,H/4-d,2*L/5,H/2]);
figure('Name',"Evolution du TEB et du TES avec SNR",'Position',[2*d,H/4-2*d,2*L/5,H/2]);
figure('Name',"Diagramme de l'oeil",'Position',[3*d,H/4-3*d,2*L/5,H/2]);
figure('Name',"Réponses impultionnelles",'Position',[4*d,H/4-4*d,2*L/5,H/2]);
figure('Name',"Réponses en fréquences",'Position',[5*d,H/4-5*d,2*L/5,H/2]);
figure('Name',"Constellation",'Position',[6*d,H/4-6*d,2*L/5,H/2]);



% Définition des constantes
Fe = 24000;
Te = 1/Fe;
Rb = 3000;
Tb = 1/Rb;
Nb = 10000;
fp = 2000;
N = 10;
K = 10;


% Simuler avec différent SNR
SNR_dBs = 0:0.1:10;
% SNR_dBs = [1000];           % pour retirer le bruit
SNRs = 10.^(SNR_dBs/10);


% Constantes des mappings
M = 2;
Ts = Tb*log2(M);
Ns = Ts/Te;


% Génération des bits
bits = [1 1 1 0 0 1 0];
bits = randi([0, 1], 1, Nb);      % pour calculer les TEB
Nb = length(bits);

x = 2 * bits - 1;


% Génération des diracs et de l'echelle temporelle
somme_dirac = kron(x, [1 zeros(1, Ns-1)]);
Temps = Te*(1:length(somme_dirac));


% Génération des diracs et de l'echelle fréquencielle
Frequences = Fe*(1/length(Temps):1/length(Temps):1);


%% Modulation

% Constantes de modulation
n0 = Ns;
n0_eg = 2*Ns;
alpha0 = 1;
alpha1 = 0.5;


% Définition des filtres
h = ones(1,Ns);
hc = [alpha0 zeros(1, length(h) - 1) alpha1];
hr = ones(1,Ns);


% Calcule de Z et C
Y0 = [1 zeros(1,N-1)];
y_dirac = kron(Y0, [1 zeros(1, Ns-1)]);
y_forme = filter(h, 1, y_dirac);
y_multi = filter(hc, 1, y_forme);
y_reception = filter(hr, 1, y_multi);

echantillonage = n0:Ns:n0+(N/log2(M)-1)*Ns;
y_echantillone = y_reception(echantillonage);

% Définition de Z
Z = zeros(K, N);
for j = 1:N
    zeros(j-1,1);
    y_echantillone(1:K-j+1);
    Z(:,j) = [zeros(1,j-1) y_echantillone(1:K-j+1)];
end

% Z = toeplitz(y_echantillone);

% Définition de C
if K == N
    C = Z\Y0';
else
    C = (Z'*Z)\Z'*Y0;
end
heg = C';


% Calcule des signaux filtés
TEBs = zeros(3,length(SNR_dBs));

for i = 1:length(SNR_dBs)
    SNR_dB = SNR_dBs(i);
    SNR = SNRs(i);
    
    
    % Filtre de mise en forme
    x_forme = filter(h, 1, somme_dirac);
    
    
    % Filtre de mise en forme
    x_multi = filter(hc, 1, x_forme);
    
    
    % Génération du bruit complexe
    Px = mean(abs(x_multi).^2);
    sigma = sqrt(Px*Ns/(2*log2(M)*SNR));
    bruit_i = sigma*randn(1,length(x_multi));
    bruit_q = sigma*randn(1,length(x_multi));
    bruit = bruit_i + 1i * bruit_q;
    
    Px = mean(abs(x_forme).^2);
    sigma = sqrt(Px*Ns/(2*log2(M)*SNR));
    bruit_i = sigma*randn(1,length(x_forme));
    bruit_q = sigma*randn(1,length(x_forme));
    bruit2 = bruit_i + 1i * bruit_q;
    
    
    % Ajout du bruit
    x_bruit = x_multi + bruit;
    x_bruit2 = x_forme + bruit2;
    
    
    % Filtre de réception
    x_reception = filter(hr, 1, x_bruit);
    x_reception2 = filter(hr, 1, x_bruit2);
    
    
    % Echantillonage
    echantillonage = n0:Ns:n0+(Nb/log2(M)-1)*Ns;
    x_echantillonne = x_reception(echantillonage);
    x_echantillonne2 = x_reception2(echantillonage);
    
    
    % Filtre d'égalisation
    x_egalisation = filter(heg, 1, x_echantillonne);
    
    
    % Reconstitution des bits
    symboles = [-1 1];
    
    [~,indices_symboles] = min(abs(x_echantillonne-symboles(:)));
    [~,indices_symboles2] = min(abs(x_echantillonne2-symboles(:)));
    [~,indices_symboles3] = min(abs(x_egalisation-symboles(:)));
    symboles_retrouves = symboles(indices_symboles);
    symboles_retrouves2 = symboles(indices_symboles2);
    symboles_retrouves3 = symboles(indices_symboles3);
    
    bits_retrouves = (symboles_retrouves + 1) / 2;
    bits_retrouves2 = (symboles_retrouves2 + 1) / 2;
    bits_retrouves3 = (symboles_retrouves3 + 1) / 2;
    
    
    % Calcule des TES et TEB
    TEB = length(find(bits_retrouves-bits ~= 0))/Nb;
    TEBs(1,i) = TEB;
    TEB2 = length(find(bits_retrouves2-bits ~= 0))/Nb;
    TEBs(2,i) = TEB2;
    TEB3 = length(find(bits_retrouves3-bits ~= 0))/Nb;
    TEBs(3,i) = TEB3;
    
end

TEB

% Courbe signaux
figure(1);
plot(Temps, real(x_forme), 'b', Temps, real(x_reception), 'r', Temps(echantillonage), real(x_egalisation), 'g');
legend('Signal formé','Signal reçu','Signal égalisé');
xlabel("Temps (en secondes)");
ylabel("Signaux");


% Courbe DSP
figure(2);
semilogy(Frequences, smoothdata(abs(fft(x_forme))));
legend('DSP du signal');
xlabel("Fréquences (en hertz)");
ylabel("DSP");


% Courbe de TEB/TES
TEB_th = 0.5 * qfunc(5/2*sqrt(10.^(SNR_dBs/10)/2)) + 0.5 * qfunc(4/5*sqrt(10.^(SNR_dBs/10)/2));
 
figure(3);
semilogy(SNR_dBs, TEBs(1,:), 'b', SNR_dBs, TEB_th, '-.b', SNR_dBs, TEBs(2,:), 'r', SNR_dBs, TEBs(3,:), 'g');
legend('TEB','TEB théorique','TEB sans canal','TEB égalisé');
xlabel("SNR (en dB)");
ylabel("TEB (en %)");


% Diagrammes de l'oeil
nb_Ns = min(Nb, 1000);

figure(4);
oeil_eg = [x_egalisation(1:end-1); x_egalisation(2:end)];
plot(Te*(0:Ns-1), real(reshape(x_reception(2*Ns+1:nb_Ns*Ns/log2(M)-Ns),Ns,[])), 'b', [0 Te*(Ns-1)], real(oeil_eg), 'r');
legend('signal receptionné','signal égalisé');
xlabel("Temps (s)");


% Filtres
figure(5);
gr = conv(h,hr);
gc = conv(gr, hc);
geg = conv(gc(1:Ns:end), heg);
plot(1:length(gr), gr, 'b', 1:length(gc), gc, 'r', 1:Ns:Ns*length(geg), geg, 'k');
legend('h * hr','h * hc * hr','h * hc * hr * heg');
xlabel("Temps (s)");
ylabel("Réponses impultionnelles");

figure(6);
hold on;
Hc = abs(fft(hc,256));
Heg = abs(fft(heg,256));
plot((1:256) /256 *Fe -Fe/2,Hc); plot((1:256) /256 *Fe -Fe/2,Heg); plot((1:256) /256 *Fe -Fe/2,Hc.*Heg);
legend('Hc','Heg','Hc * Heg');
xlabel("Fréquences (Hz)");
ylabel("Réponses en fréquences");


% Constellation
figure(7);
plot(real(x_echantillonne), imag(x_echantillonne), '*b', real(x_egalisation), imag(x_egalisation), '*g', [-1 1], [0 0], '*r');
legend('constellation de réception','contellation égalisation','constellation initiale');
