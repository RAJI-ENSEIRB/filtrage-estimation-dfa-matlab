function [P_moy, f] = periodogram_welch(signal, fs, L, overlap)
    % signal : signal à analyser
    % fs : fréquence d'échantillonnage
    % L : Longueur des segments
    % overlap : chevauchement (en nombre d'échantillons)

    N = length(signal);
    step = L - overlap; % Pas entre les segments
    K = floor((N - overlap) / step); % Nombre de segments
    P_moy = zeros(L/2, 1); % Initialiser le périodogramme moyen
    window = hanning(L); % Fenêtre de Hanning

    % Normalisation de la fenêtre pour conserver la puissance
    U = sum(window.^2) / L; % Facteur de normalisation

    for k = 1:K
        start_idx = (k-1)*step + 1;
        end_idx = start_idx + L - 1;
        if end_idx > N
            break; % Éviter de dépasser la longueur du signal
        end
        segment = signal(start_idx:end_idx); % Extraire un segment chevauchant
        segment = segment .* window; % Appliquer la fenêtre de Hanning
        [P_seg, ~] = periodogram_brut(segment, fs); % Calculer le périodogramme du segment
        P_moy = P_moy + P_seg; % Accumuler le périodogramme
    end

    P_moy = P_moy / K / U; % Moyenne sur tous les segments et normalisation
    f = (0:(L/2-1))*(fs/L); % Fréquences correspondantes
end
