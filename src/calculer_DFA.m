function [H, alpha, F2, N_values] = calculer_DFA(signal, degree)
    % Normalisation du signal
    signal = (signal - mean(signal)) / std(signal);
    M = length(signal);
    
    % Tailles de fenêtres
    N_min = max(6, degree + 2);
    N_max = floor(M/4);
    N_values = unique(round(logspace(log10(N_min), log10(N_max), 20)));
    
    % Profil intégré
    profil = cumsum(signal - mean(signal));
    
    % Initialisation
    F2 = zeros(size(N_values));
    
    % Pour chaque taille de fenêtre
    for i = 1:length(N_values)
        N = N_values(i);
        L = floor(M/N);
        F2_temp = 0;
        
        for j = 1:L
            debut = (j-1)*N + 1;
            fin = j*N;
            segment = profil(debut:fin);
            
            % Régression polynomiale
            t = (1:N)';
            if degree == 1
                X = [ones(N,1), t];
            else  % degree == 2
                X = [ones(N,1), t, t.^2];
            end
            
            % Calcul de la tendance
            beta = X\segment;
            tendance = X*beta;
            
            % Accumulation de F2
            F2_temp = F2_temp + sum((segment - tendance).^2);
        end
        
        F2(i) = sqrt(F2_temp/(L*N));
    end
    
    % Régression log-log
    logN = log(N_values);
    logF2 = log(F2);
    X = [ones(length(N_values),1), logN'];
    coeffs = X\logF2';
    alpha = coeffs(2);
    H = 1 - alpha;
    
    % Affichage
    figure;
    subplot(2,1,1);
    plot(signal);
    title(sprintf('Signal original - DFA degré %d', degree));
    xlabel('Échantillons');
    ylabel('Amplitude');
    grid on;
    
    subplot(2,1,2);
    loglog(N_values, F2, 'b.', 'MarkerSize', 15);
    hold on;
    loglog(N_values, exp(coeffs(1) + coeffs(2)*logN), 'r-', 'LineWidth', 1.5);
    title(sprintf('DFA: H = %.3f, α = %.3f (degré %d)', H, alpha, degree));
    xlabel('log(N)');
    ylabel('log(F^2(N))');
    grid on;
    legend('Données', 'Régression', 'Location', 'best');
end