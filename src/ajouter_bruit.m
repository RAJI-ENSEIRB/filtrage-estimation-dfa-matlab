function signal_bruite = ajouter_bruit(signal, RSB_dB)
    % Fonction pour ajouter du bruit blanc Gaussien à un signal pour obtenir un RSB donné
    % signal : le signal d'entrée (vecteur)
    % RSB_dB : rapport signal/bruit en dB (ex. 5, 10, 15 dB)

    % 1. Calcul de la puissance du signal (moyenne des carrés)
    P_signal = mean(signal.^2);

    % 2. Calcul de la puissance du bruit pour obtenir le RSB souhaité
    RSB = 10^(RSB_dB / 10); % Conversion de dB en ratio linéaire
    P_bruit = P_signal / RSB; % Puissance du bruit à ajouter

    % 3. Génération du bruit blanc Gaussien
    bruit = sqrt(P_bruit) * randn(size(signal));

    % 4. Ajout du bruit au signal
    signal_bruite = signal + bruit;
end
