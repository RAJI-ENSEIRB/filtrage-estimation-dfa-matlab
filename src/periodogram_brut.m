function [P, f] = periodogram_brut(signal, fs)
    N = length(signal); % Longueur du signal
    X = fft(signal); % Transformée de Fourier
    P = (1/N) * abs(X).^2; % Périodogramme brut (densité spectrale de puissance)
    f = (0:N-1)*(fs/N); % Fréquences correspondantes
    P = P(1:floor(N/2)); % On ne garde que les fréquences positives
    f = f(1:floor(N/2)); % On ne garde que les fréquences positives
end
