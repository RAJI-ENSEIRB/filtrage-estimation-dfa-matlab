function [P_moy, f] = periodogram_bartlett(signal, fs, L)
    % signal : signal à analyser
    % fs : fréquence d'échantillonnage
    % L : Longueur des segments

    N = length(signal);
    K = floor(N / L); % Nombre de segments
    P_moy = zeros(L/2, 1); % Initialiser le périodogramme moyen

    for k = 1:K
        segment = signal((k-1)*L + (1:L)); % Extraire un segment
        [P_seg, f] = periodogram_brut(segment, fs); % Calculer le périodogramme du segment
        P_moy = P_moy + P_seg; % Accumuler le périodogramme
    end
    
    P_moy = P_moy / K; % Moyenne sur tous les segments
end
