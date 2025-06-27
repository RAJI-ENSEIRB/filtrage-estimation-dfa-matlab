% Script MATLAB pour la génération et l'analyse d'un bruit blanc

% Paramètres
N = 1000; % Taille du signal
sigma = 1; % Variance du bruit

% 1. Génération d’un bruit blanc de moyenne nulle et de variance sigma^2
bruit_blanc = sigma * randn(N, 1); % Génération du bruit blanc de variance sigma^2

% 2. Autocorrélation théorique, biaisée et non biaisée
% Autocorrélation théorique : nul partout sauf à k = 0, où elle vaut sigma^2
autocorr_theorique = sigma^2 * [1, zeros(1, N-1)];

% Estimation de l'autocorrélation biaisée
[autocorr_biaisee, lags] = xcorr(bruit_blanc, 'biased');

% Estimation de l'autocorrélation non biaisée
[autocorr_non_biaisee, ~] = xcorr(bruit_blanc, 'unbiased');

% Tracé de l'autocorrélation théorique, biaisée et non biaisée
figure;
plot(lags, autocorr_biaisee, 'b', 'DisplayName', 'Autocorrélation biaisée');
hold on;
plot(lags, autocorr_non_biaisee, 'g', 'DisplayName', 'Autocorrélation non biaisée');
stem(0, autocorr_theorique(1), 'r', 'DisplayName', 'Autocorrélation théorique');
title('Autocorrélation du bruit blanc');
legend show;
xlabel('Décalage');
ylabel('Autocorrélation');
hold off;

% 3. Spectre de puissance d’une réalisation de ce bruit blanc
% Calcul du spectre de puissance avec la transformée de Fourier
spectre_bruit = abs(fft(bruit_blanc)).^2 / N;
freqs = linspace(0, 1, N); % Fréquences normalisées

% Tracé du spectre de puissance
figure;
plot(freqs, spectre_bruit);
title('Spectre de puissance du bruit blanc');
xlabel('Fréquence normalisée');
ylabel('Puissance');

% 4. Densité spectrale de puissance (DSP) estimée avec l’algorithme de Welch
[pxx, f] = pwelch(bruit_blanc, [], [], [], 1); % DSP avec Welch

% Tracé de la DSP
figure;
plot(f, 10*log10(pxx));
title('Densité spectrale de puissance (DSP) du bruit blanc');
xlabel('Fréquence normalisée');
ylabel('Puissance (dB)');

%% Analyse et Programmation des Périodogrammes de Daniell, Bartlett, et Welch

% Paramètres du signal
fs = 1000;    % Fréquence d'échantillonnage (Hz)
N = 1000;     % Taille du signal
sigma = 1;    % Variance du bruit blanc

% Génération du bruit blanc Gaussien centré
bruit_blanc = sigma * randn(N, 1);

% Calcul du périodogramme brut (DSP empirique)
[P_brut, f_brut] = periodogram_brut(bruit_blanc, fs);

% Calcul du périodogramme de Daniell (avec une fenêtre de lissage M = 5)
P_daniell = periodogram_daniell(P_brut, 5);

% Calcul du périodogramme de Bartlett (segments non chevauchants, longueur 128)
[P_bartlett, f_bartlett] = periodogram_bartlett(bruit_blanc, fs, 128);

% Calcul du périodogramme de Welch (segments chevauchants, longueur 128, chevauchement 64)
[P_welch, f_welch] = periodogram_welch(bruit_blanc, fs, 128, 64);

% Tracé de la DSP théorique (constante à sigma^2)
DSP_theorique = sigma^2 * ones(size(f_brut)); % Valeur constante égale à la variance du bruit blanc

% ---- Affichage de toutes les méthodes sur un même graphe ----
figure;
plot(f_brut, 10*log10(P_brut), 'b', 'DisplayName', 'Périodogramme brut');
hold on;
plot(f_brut, 10*log10(P_daniell), 'r', 'DisplayName', 'Périodogramme de Daniell');
plot(f_bartlett, 10*log10(P_bartlett), 'g', 'DisplayName', 'Périodogramme de Bartlett');
plot(f_welch, 10*log10(P_welch), 'm', 'DisplayName', 'Périodogramme de Welch');
plot(f_brut, 10*log10(DSP_theorique), 'k--', 'DisplayName', 'DSP théorique');
title('Comparaison des Périodogrammes (Brut, Daniell, Bartlett, Welch) et DSP théorique');
xlabel('Fréquence (Hz)');
ylabel('Puissance (dB)');
legend;
hold off;

% ---- Tracé de chaque périodogramme individuellement avec la DSP théorique ----

% Périodogramme brut et DSP théorique
figure;
plot(f_brut, 10*log10(P_brut), 'b', 'DisplayName', 'Périodogramme brut');
hold on;
plot(f_brut, 10*log10(DSP_theorique), 'k--', 'DisplayName', 'DSP théorique');
title('Périodogramme brut et DSP théorique');
xlabel('Fréquence (Hz)');
ylabel('Puissance (dB)');
legend;
hold off;

% Périodogramme de Daniell et DSP théorique
figure;
plot(f_brut, 10*log10(P_daniell), 'r', 'DisplayName', 'Périodogramme de Daniell');
hold on;
plot(f_brut, 10*log10(DSP_theorique), 'k--', 'DisplayName', 'DSP théorique');
title('Périodogramme de Daniell et DSP théorique');
xlabel('Fréquence (Hz)');
ylabel('Puissance (dB)');
legend;
hold off;

% Périodogramme de Bartlett et DSP théorique
figure;
plot(f_bartlett, 10*log10(P_bartlett), 'g', 'DisplayName', 'Périodogramme de Bartlett');
hold on;
plot(f_brut, 10*log10(DSP_theorique), 'k--', 'DisplayName', 'DSP théorique');
title('Périodogramme de Bartlett et DSP théorique');
xlabel('Fréquence (Hz)');
ylabel('Puissance (dB)');
legend;
hold off;

% Périodogramme de Welch et DSP théorique
figure;
plot(f_welch, 10*log10(P_welch), 'm', 'DisplayName', 'Périodogramme de Welch');
hold on;
plot(f_brut, 10*log10(DSP_theorique), 'k--', 'DisplayName', 'DSP théorique');
title('Périodogramme de Welch et DSP théorique');
xlabel('Fréquence (Hz)');
ylabel('Puissance (dB)');
legend;
hold off;

%% Étude et Programmation du Corrélogramme
% Paramètres
fs = 1000;   % Fréquence d'échantillonnage (Hz)
N = 1000;    % Longueur du signal
sigma = 1;   % Variance du bruit blanc

% Génération du bruit blanc Gaussien centré
bruit_blanc = sigma * randn(N, 1);

% 1. Calcul du corrélogramme
% Estimation de l'autocorrélation
Rxx = xcorr(bruit_blanc, 'biased'); % Autocorrélation normalisée
Rxx = Rxx(N:end); % Garder uniquement les lags positifs

% Calcul de la densité spectrale de puissance via la transformée de Fourier
P_corr = abs(fft(Rxx)); % Magnitude de la FFT de l'autocorrélation
P_corr = P_corr(1:N/2); % Prendre seulement les fréquences positives
f_corr = (0:N/2-1)*(fs/N); % Fréquences associées

% 2. Calcul du périodogramme brut avec la fonction periodogram_brut
[P_brut, f_brut] = periodogram_brut(bruit_blanc, fs);

% 3. DSP théorique (constante à sigma^2 pour le bruit blanc)
DSP_theorique = sigma^2 * ones(size(f_brut));

% Comparaison des résultats
figure;
plot(f_corr, 10*log10(P_corr), 'b', 'DisplayName', 'Corrélogramme');
hold on;
plot(f_brut, 10*log10(P_brut), 'r', 'DisplayName', 'Périodogramme brut');
plot(f_brut, 10*log10(DSP_theorique), 'k--', 'DisplayName', 'DSP théorique');
title('Comparaison Corrélogramme, Périodogramme brut et DSP théorique');
xlabel('Fréquence (Hz)');
ylabel('Puissance (dB)');
legend;
hold off;


%% Ajout de Bruit Additif à un Signal avec un RSB Donné : Implémentation pour des Signaux de Weierstrass et de Parole sous MATLAB

% Chargement du signal de Weierstrass depuis le fichier fourni
data = load('data_Weierstrass.mat'); % Charger le fichier .mat
signal_weierstrass = data.data{1, 1}; 

% Vérifiez si le signal est un vecteur (sinon, le transformer)
if iscell(signal_weierstrass)
    signal_weierstrass = signal_weierstrass{1};
end

% Vérifiez les dimensions du signal
signal_weierstrass = signal_weierstrass(:); % Assurez-vous qu'il s'agit d'un vecteur colonne

% Définir un vecteur temporel x basé sur la taille du signal
N = length(signal_weierstrass);
x = linspace(0, 1, N); % Temps normalisé entre 0 et 1

% Tester avec différents RSB
RSB_dB = [5, 10, 15]; % RSB à tester en dB

% Bruitage du signal de Weierstrass pour chaque RSB
for i = 1:length(RSB_dB)
    signal_bruite = ajouter_bruit(signal_weierstrass, RSB_dB(i));
    
    % Affichage des résultats
    figure;
    plot(x, signal_weierstrass, 'b', 'DisplayName', 'Signal original');
    hold on;
    plot(x, signal_bruite, 'r', 'DisplayName', ['Signal bruité (RSB = ', num2str(RSB_dB(i)), ' dB)']);
    title(['Signal de Weierstrass avec RSB = ', num2str(RSB_dB(i)), ' dB']);
    xlabel('Temps');
    ylabel('Amplitude');
    legend show;
    hold off;
end


% Chargement du fichier .mat
data = load('fcno03fz.mat'); % Charger le fichier .mat

% Extraction du signal à partir de la variable 'fcno03fz'
signal_parole = data.fcno03fz;

% Définir une fréquence d'échantillonnage si nécessaire (ex. : 16 kHz)
fs = 16000; % Supposition d'une fréquence d'échantillonnage

% Tester avec différents RSB
RSB_dB = [5, 10, 15]; % RSB à tester en dB

for i = 1:length(RSB_dB)
    signal_bruite_parole = ajouter_bruit(signal_parole, RSB_dB(i));
    
    % Écouter le signal bruité
    sound(signal_bruite_parole, fs); % Jouer le son avec la fréquence d'échantillonnage fs
    pause(5); % Pause pour laisser le temps d'écouter

    % Affichage de l'onde
    t = (0:length(signal_parole)-1) / fs; % Temps
    figure;
    plot(t, signal_parole, 'b', 'DisplayName', 'Signal original');
    hold on;
    plot(t, signal_bruite_parole, 'r', 'DisplayName', ['Signal bruité (RSB = ', num2str(RSB_dB(i)), ' dB)']);
    title(['Signal de parole avec RSB = ', num2str(RSB_dB(i)), ' dB']);
    xlabel('Temps (s)');
    ylabel('Amplitude');
    legend show;
    hold off;
end
%% Analyse Temporelle et Spectrale d'un Signal : Comparaison entre Signal Original et Signal Bruité
% Chargement du signal Weierstrass
data = load('data_Weierstrass.mat');
signal_weierstrass = data.data{1, 1};
if iscell(signal_weierstrass)
    signal_weierstrass = signal_weierstrass{1};
end
signal_weierstrass = signal_weierstrass(:);

% Paramètres
window = hamming(256);
noverlap = 128;
nfft = 512;
fs = 1; % Fréquence d'échantillonnage normalisée
t = (0:length(signal_weierstrass)-1)/fs;

% Pour chaque RSB
RSB_dB = [5, 10, 15];
for i = 1:length(RSB_dB)
    signal_bruite = ajouter_bruit(signal_weierstrass, RSB_dB(i));
    
    % Figure 1: Signal original
    figure('Position', [100 100 800 600]);
    ax1 = subplot(2,1,1);
    plot(t, signal_weierstrass, 'LineWidth', 1);
    title('Signal original');
    ylabel('Amplitude');
    xlabel('Temps (s)');
    grid on;
    
    ax2 = subplot(2,1,2);
    [s,f,t_spec] = spectrogram(signal_weierstrass, window, noverlap, nfft, fs, 'yaxis');
    surf(t_spec, f, 10*log10(abs(s)), 'EdgeColor', 'none');
    axis tight;
    view(0,90);
    title('Spectrogramme du signal original');
    ylabel('Fréquence normalisée');
    xlabel('Temps (s)');
    colorbar;
    linkaxes([ax1 ax2], 'x');
    
    % Figure 2: Signal bruité
    figure('Position', [100 100 800 600]);
    ax3 = subplot(2,1,1);
    plot(t, signal_bruite, 'LineWidth', 1);
    title(['Signal bruité (RSB = ' num2str(RSB_dB(i)) ' dB)']);
    ylabel('Amplitude');
    xlabel('Temps (s)');
    grid on;
    
    ax4 = subplot(2,1,2);
    [s,f,t_spec] = spectrogram(signal_bruite, window, noverlap, nfft, fs, 'yaxis');
    surf(t_spec, f, 10*log10(abs(s)), 'EdgeColor', 'none');
    axis tight;
    view(0,90);
    title(['Spectrogramme du signal bruité (RSB = ' num2str(RSB_dB(i)) ' dB)']);
    ylabel('Fréquence normalisée');
    xlabel('Temps (s)');
    colorbar;
    linkaxes([ax3 ax4], 'x');
end

%%
% PARTIE 2
%% les étapes de la méthode DFA

% Nettoyage
clear all;
close all;
clc;

% Chargement des données
data = load('data_Weierstrass.mat');
nb_regularites = size(data.data, 1);

% Initialisation
H_values_deg1 = zeros(nb_regularites, 1);
H_values_deg2 = zeros(nb_regularites, 1);
regularites_theoriques = (1:nb_regularites)*0.1;

% Analyse pour chaque régularité
for i = 1:nb_regularites
    signal = data.data{i,1};
    
    % Si c'est la première régularité, on fait les visualisations détaillées
    if i == 1
        % 1. Signal initial et centré
        figure('Position', [100 100 1200 800]);
        subplot(2,2,1);
        plot(signal);
        title('Signal initial');
        xlabel('Échantillons'); ylabel('Amplitude');
        grid on;
        
        signal_centre = signal - mean(signal);
        subplot(2,2,2);
        plot(signal_centre);
        title('Signal centré');
        xlabel('Échantillons'); ylabel('Amplitude');
        grid on;
        
        % 2. Profil (signal intégré)
        profil = cumsum(signal_centre);
        subplot(2,2,3);
        plot(profil);
        title('Profil (signal intégré)');
        xlabel('Échantillons'); ylabel('Amplitude');
        grid on;
        
        % Calcul des paramètres pour le découpage
        N = floor(length(signal)/10); % Taille des segments
        L = floor(length(signal)/N);  % Nombre de segments complets
        
        % 3. Exemple d'un segment avec tendances
        segment_num = floor(L/2); % On prend un segment au milieu
        debut = (segment_num-1)*N + 1;
        fin = segment_num*N;
        
        % Extraction du segment
        t = (1:N)';
        segment = profil(debut:fin);
        
        % Calcul des tendances
        % Degré 1
        X1 = [ones(N,1), t];
        beta1 = X1\segment;
        tendance1 = X1*beta1;
        
        % Degré 2
        X2 = [ones(N,1), t, t.^2];
        beta2 = X2\segment;
        tendance2 = X2*beta2;
        
        subplot(2,2,4);
        plot(debut:fin, segment, 'b.', 'MarkerSize', 10);
        hold on;
        plot(debut:fin, tendance1, 'r-', 'LineWidth', 2);
        plot(debut:fin, tendance2, 'g-', 'LineWidth', 2);
        title('Segment exemple avec tendances');
        xlabel('Échantillons'); ylabel('Amplitude');
        legend('Profil', 'Tendance deg. 1', 'Tendance deg. 2');
        grid on;
        
        % 4. Tendance globale et résidu
        figure('Position', [100 100 1200 400]);
        subplot(1,2,1);
        tendance_globale1 = zeros(L*N, 1);
        tendance_globale2 = zeros(L*N, 1);
        
        for j = 1:L
            debut = (j-1)*N + 1;
            fin = j*N;
            segment = profil(debut:fin);
            t = (1:N)';
            
            % Degré 1
            X1 = [ones(N,1), t];
            beta1 = X1\segment;
            tendance_globale1(debut:fin) = X1*beta1;
            
            % Degré 2
            X2 = [ones(N,1), t, t.^2];
            beta2 = X2\segment;
            tendance_globale2(debut:fin) = X2*beta2;
        end
        
        plot(1:L*N, profil(1:L*N));
        hold on;
        plot(1:L*N, tendance_globale1, 'r');
        plot(1:L*N, tendance_globale2, 'g');
        title('Profil et tendances globales');
        xlabel('Échantillons'); ylabel('Amplitude');
        legend('Profil', 'Tendance deg. 1', 'Tendance deg. 2');
        grid on;
        
        % Résidus
        subplot(1,2,2);
        residu1 = profil(1:L*N) - tendance_globale1;
        residu2 = profil(1:L*N) - tendance_globale2;
        plot(1:L*N, residu1, 'r');
        hold on;
        plot(1:L*N, residu2, 'g');
        title('Résidus (Profil - Tendance)');
        xlabel('Échantillons'); ylabel('Amplitude');
        legend('Résidu deg. 1', 'Résidu deg. 2');
        grid on;
        
        % 5. Nouvelle figure : Analyse log-log de F²(N)
        % Calcul de F²(N) pour différentes tailles de fenêtre
        N_min = 10;
        N_max = floor(length(signal)/4);
        N_values = unique(round(logspace(log10(N_min), log10(N_max), 20)));
        
        F2_deg1 = zeros(size(N_values));
        F2_deg2 = zeros(size(N_values));
        
        for k = 1:length(N_values)
            N_curr = N_values(k);
            L_curr = floor(length(signal)/N_curr);
            F2_temp1 = 0;
            F2_temp2 = 0;
            
            for j = 1:L_curr
                debut = (j-1)*N_curr + 1;
                fin = j*N_curr;
                segment = profil(debut:fin);
                t = (1:N_curr)';
                
                % Degré 1
                X1 = [ones(N_curr,1), t];
                beta1 = X1\segment;
                tendance1 = X1*beta1;
                F2_temp1 = F2_temp1 + sum((segment - tendance1).^2);
                
                % Degré 2
                X2 = [ones(N_curr,1), t, t.^2];
                beta2 = X2\segment;
                tendance2 = X2*beta2;
                F2_temp2 = F2_temp2 + sum((segment - tendance2).^2);
            end
            
            F2_deg1(k) = sqrt(F2_temp1/(L_curr*N_curr));
            F2_deg2(k) = sqrt(F2_temp2/(L_curr*N_curr));
        end
        
        % Régression log-log
        logN = log(N_values);
        logF2_1 = log(F2_deg1);
        logF2_2 = log(F2_deg2);
        
        % Coefficients de régression
        X_reg = [ones(length(N_values),1), logN'];
        coeffs1 = X_reg\logF2_1';
        coeffs2 = X_reg\logF2_2';
        
        alpha1 = coeffs1(2);
        alpha2 = coeffs2(2);
        H1 = alpha1 - 1;
        H2 = alpha2 - 1;
        
        % Affichage de la régression log-log
        figure('Position', [100 600 800 400]);
        loglog(N_values, F2_deg1, 'ro', 'MarkerSize', 8);
        hold on;
        loglog(N_values, F2_deg2, 'go', 'MarkerSize', 8);
        loglog(N_values, exp(coeffs1(1) + coeffs1(2)*logN), 'r--', 'LineWidth', 1.5);
        loglog(N_values, exp(coeffs2(1) + coeffs2(2)*logN), 'g--', 'LineWidth', 1.5);
        title(sprintf('Analyse log-log F²(N) vs N\nH1=%.3f, α1=%.3f, H2=%.3f, α2=%.3f', H1, alpha1, H2, alpha2));
        xlabel('log(N)');
        ylabel('log(F²(N))');
        legend('Données deg. 1', 'Données deg. 2', 'Régression deg. 1', 'Régression deg. 2', 'Location', 'best');
        grid on;
    end
    
    % DFA avec degré 1
    [H1, alpha1] = calculer_DFA(signal, 1);
    H_values_deg1(i) = H1;
    
    % DFA avec degré 2
    [H2, alpha2] = calculer_DFA(signal, 2);
    H_values_deg2(i) = H2;
    
    fprintf('Régularité %.1f :\n', regularites_theoriques(i));
    fprintf('  Degré 1 : H = %.3f, alpha = %.3f\n', H1, alpha1);
    fprintf('  Degré 2 : H = %.3f, alpha = %.3f\n\n', H2, alpha2);
end

% Comparaison des résultats
figure('Position', [100 100 600 400]);
plot(regularites_theoriques, H_values_deg1, 'bo-', 'LineWidth', 1.5);
hold on;
plot(regularites_theoriques, H_values_deg2, 'ro-', 'LineWidth', 1.5);
plot(regularites_theoriques, regularites_theoriques, 'k--');
xlabel('Régularité théorique');
ylabel('H estimé');
title('Comparaison des estimations DFA selon le degré polynomial');
legend('Degré 1', 'Degré 2', 'Valeur théorique', 'Location', 'best');
grid on;

% Tableau récapitulatif
resultats = table(regularites_theoriques', H_values_deg1, H_values_deg2, ...
    'VariableNames', {'Regularite', 'H_deg1', 'H_deg2'});
disp('Résultats par régularité :');
disp(resultats);
%%  une réalisation d’un bruit blanc Gaussien et une réalisation d'un bruit rose
% Nettoyage
clear all;
close all;
clc;

% Paramètres
N = 10000;           % Nombre d'échantillons
t = (1:N)';         % Vecteur temps

% Création du bruit blanc gaussien
bruit_blanc = randn(N, 1);

% Création du bruit rose (méthode améliorée)
nfft = N;
f = (0:nfft/2)';
P_rose = 1./(f + eps).^1;  % Spectre en 1/f^1
P_rose(1) = P_rose(2);     % Éviter la singularité à f=0
P_rose = P_rose/sqrt(sum(P_rose.^2)); % Normalisation
P = [P_rose; flipud(P_rose(2:end-1))];

% Génération des phases aléatoires
phases = exp(2i*pi*rand(nfft, 1));
phases(1) = 1;
phases(nfft/2+1) = 1;
phases(nfft/2+2:end) = conj(flipud(phases(2:nfft/2)));

% Création du signal dans le domaine fréquentiel
X = sqrt(P) .* phases;
bruit_rose = real(ifft(X));

% Normalisation des signaux
bruit_blanc = (bruit_blanc - mean(bruit_blanc)) / std(bruit_blanc);
bruit_rose = (bruit_rose - mean(bruit_rose)) / std(bruit_rose);

% Définition des tendances
tendance_exp = 2 * exp(t/N * 3);          % Tendance exponentielle
tendance_log = 3 * log(t + 1);            % Tendance logarithmique

% Ajout des tendances aux bruits
bruit_blanc_exp = bruit_blanc + tendance_exp;
bruit_blanc_log = bruit_blanc + tendance_log;
bruit_rose_exp = bruit_rose + tendance_exp;
bruit_rose_log = bruit_rose + tendance_log;

% Visualisation des signaux
figure('Position', [100 100 1200 800]);

% Bruit blanc et ses variantes
subplot(3,2,1);
plot(bruit_blanc);
title('Bruit blanc pur');
xlabel('Échantillons');
ylabel('Amplitude');
grid on;

subplot(3,2,2);
plot(bruit_blanc_exp);
title('Bruit blanc + tendance exp.');
xlabel('Échantillons');
ylabel('Amplitude');
grid on;

subplot(3,2,3);
plot(bruit_rose);
title('Bruit rose pur');
xlabel('Échantillons');
ylabel('Amplitude');
grid on;

subplot(3,2,4);
plot(bruit_rose_exp);
title('Bruit rose + tendance exp.');
xlabel('Échantillons');
ylabel('Amplitude');
grid on;

subplot(3,2,5);
plot(bruit_blanc_log);
title('Bruit blanc + tendance log.');
xlabel('Échantillons');
ylabel('Amplitude');
grid on;

subplot(3,2,6);
plot(bruit_rose_log);
title('Bruit rose + tendance log.');
xlabel('Échantillons');
ylabel('Amplitude');
grid on;

% Application du DFA
fprintf('=== Résultats DFA ===\n\n');

fprintf('Bruit blanc pur:\n');
[H_bw, alpha_bw] = calculer_DFA(bruit_blanc, 1);
fprintf('H = %.3f, alpha = %.3f\n\n', H_bw, alpha_bw);

fprintf('Bruit rose pur:\n');
[H_br, alpha_br] = calculer_DFA(bruit_rose, 1);
fprintf('H = %.3f, alpha = %.3f\n\n', H_br, alpha_br);

fprintf('Bruit blanc + tendance exponentielle:\n');
[H_bw_exp, alpha_bw_exp] = calculer_DFA(bruit_blanc_exp, 1);
fprintf('H = %.3f, alpha = %.3f\n\n', H_bw_exp, alpha_bw_exp);

fprintf('Bruit blanc + tendance logarithmique:\n');
[H_bw_log, alpha_bw_log] = calculer_DFA(bruit_blanc_log, 1);
fprintf('H = %.3f, alpha = %.3f\n\n', H_bw_log, alpha_bw_log);

fprintf('Bruit rose + tendance exponentielle:\n');
[H_br_exp, alpha_br_exp] = calculer_DFA(bruit_rose_exp, 1);
fprintf('H = %.3f, alpha = %.3f\n\n', H_br_exp, alpha_br_exp);

fprintf('Bruit rose + tendance logarithmique:\n');
[H_br_log, alpha_br_log] = calculer_DFA(bruit_rose_log, 1);
fprintf('H = %.3f, alpha = %.3f\n\n', H_br_log, alpha_br_log);

% Création d'un tableau récapitulatif
Types = {'Bruit blanc pur'; 'Bruit rose pur'; 
         'Bruit blanc + exp'; 'Bruit blanc + log';
         'Bruit rose + exp'; 'Bruit rose + log'};
H_values = [H_bw; H_br; H_bw_exp; H_bw_log; H_br_exp; H_br_log];
Alpha_values = [alpha_bw; alpha_br; alpha_bw_exp; alpha_bw_log; alpha_br_exp; alpha_br_log];

resultats = table(Types, H_values, Alpha_values, ...
    'VariableNames', {'Type_Signal', 'H', 'Alpha'});
disp('Tableau récapitulatif :');
disp(resultats);

% Sauvegarde des résultats dans un fichier
save('resultats_DFA.mat', 'resultats', 'H_values', 'Alpha_values');

%% la base de données d’enregistrements vocaux RAVDESS
% Nettoyage
clear all;
close all;
clc;

% Dossier principal
dossier_principal = 'audio_speech_actors_01-24';

% 1. Analyse par émotion
acteur = '01';    % Acteur choisi
phrase = '01';    % Phrase choisie
emotions = {'01', '03', '04', '05', '06'};  % Neutre, Heureux, Triste, Colère, Peur
emotions_names = {'Neutre', 'Heureux', 'Triste', 'Colère', 'Peur'};

% Initialisation pour les émotions
H_emotions = zeros(length(emotions), 1);
alpha_emotions = zeros(length(emotions), 1);
found_files = false(length(emotions), 1);

% Analyse pour chaque émotion
for i = 1:length(emotions)
    % Construction du chemin du fichier
    filename = ['03-01-' emotions{i} '-01-' phrase '-01-' acteur '.wav'];
    filepath = fullfile(dossier_principal, ['Actor_' acteur], filename);
    
    % Vérification de l'existence du fichier
    if ~exist(filepath, 'file')
        warning(['Fichier non trouvé: ' filepath]);
        continue;
    end
    
    found_files(i) = true;
    [y, Fs] = audioread(filepath);
    signal = y(:,1);
    
    % Calcul DFA pour l'émotion
    [H, alpha] = calculer_DFA(signal, 1);
    H_emotions(i) = H;
    alpha_emotions(i) = alpha;
end

% 2. Analyse par genre
acteurs_homme = {'01', '03', '05', '07', '09'};
acteurs_femme = {'02', '04', '06', '08', '10'};
emotion = '01';  % Émotion neutre pour la comparaison des genres

% Initialisation pour les genres
H_hommes = zeros(length(acteurs_homme), 1);
H_femmes = zeros(length(acteurs_femme), 1);
found_hommes = false(length(acteurs_homme), 1);
found_femmes = false(length(acteurs_femme), 1);

% Figure pour la comparaison des genres
figure('Position', [100 100 800 400]);

% Analyse des voix masculines
for i = 1:length(acteurs_homme)
    filename = ['03-01-' emotion '-01-' phrase '-01-' acteurs_homme{i} '.wav'];
    filepath = fullfile(dossier_principal, ['Actor_' acteurs_homme{i}], filename);
    
    if ~exist(filepath, 'file')
        warning(['Fichier non trouvé: ' filepath]);
        continue;
    end
    
    found_hommes(i) = true;
    [y, Fs] = audioread(filepath);
    signal = y(:,1);
    
    [H, alpha] = calculer_DFA(signal, 1);
    H_hommes(i) = H;
end

% Analyse des voix féminines
for i = 1:length(acteurs_femme)
    filename = ['03-01-' emotion '-01-' phrase '-01-' acteurs_femme{i} '.wav'];
    filepath = fullfile(dossier_principal, ['Actor_' acteurs_femme{i}], filename);
    
    if ~exist(filepath, 'file')
        warning(['Fichier non trouvé: ' filepath]);
        continue;
    end
    
    found_femmes(i) = true;
    [y, Fs] = audioread(filepath);
    signal = y(:,1);
    
    [H, alpha] = calculer_DFA(signal, 1);
    H_femmes(i) = H;
end

% Affichage du graphique de comparaison des genres
plot(1:length(acteurs_homme), H_hommes, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
plot(1:length(acteurs_femme), H_femmes, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
title('Comparaison des valeurs de H par genre');
xlabel('Numéro d''acteur');
ylabel('Valeur de H');
legend('Hommes', 'Femmes', 'Location', 'best');
ylim([0 1]);
xlim([0 6]);

% Figure pour comparer les émotions
figure('Position', [100 500 600 400]);
emotions_idx = 1:length(emotions);
bar(emotions_idx, H_emotions);
set(gca, 'XTick', emotions_idx);
set(gca, 'XTickLabel', emotions_names);
title('Comparaison des valeurs de H par émotion');
ylabel('Valeur de H');
grid on;
ylim([0 1]);

% Affichage des résultats par émotion
fprintf('\n=== Résultats par émotion (même acteur, même phrase) ===\n');
for i = 1:length(emotions)
    if found_files(i)
        fprintf('%s: H = %.3f, α = %.3f\n', emotions_names{i}, H_emotions(i), alpha_emotions(i));
    else
        fprintf('%s: Données non disponibles\n', emotions_names{i});
    end
end

% Affichage des résultats par genre
fprintf('\n=== Résultats par genre (même phrase, même émotion) ===\n');
if any(found_hommes)
    H_hommes_valid = H_hommes(found_hommes);
    fprintf('Hommes (moyenne): H = %.3f ± %.3f\n', mean(H_hommes_valid), std(H_hommes_valid));
else
    fprintf('Hommes: Données non disponibles\n');
end

if any(found_femmes)
    H_femmes_valid = H_femmes(found_femmes);
    fprintf('Femmes (moyenne): H = %.3f ± %.3f\n', mean(H_femmes_valid), std(H_femmes_valid));
else
    fprintf('Femmes: Données non disponibles\n');
end
% Fonction auxiliaire
function [emotion, acteur] = extraire_metadata(filename)
    parts = split(filename, '-');
    emotion = str2double(parts{3});
    acteur = str2double(parts{7}(1:end-4));
end




