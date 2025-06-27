function [P_lisse] = periodogram_daniell(P_brut, M)
    % P_brut : périodogramme brut (DSP)
    % M : Taille de la fenêtre de lissage (moyenne mobile)
    window = ones(1, M) / M; % Fenêtre de lissage
    P_lisse = conv(P_brut, window, 'same'); % Convolution pour lisser
end
