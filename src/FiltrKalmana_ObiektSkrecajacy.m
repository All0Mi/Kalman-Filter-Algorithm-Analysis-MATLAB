std_r = 10;                     % odchylenie pomiaru (z dokładnością do 10m)
std_u = 1;                      % odchylenie procesu (m/s)
dt = 1;                         % krok czasowy (s)
x = [0; 0];                     % Początkowa estymata stanu [pozycja; prędkość] (m, m/s)
P = eye(2);                     % Początkowa estymata kowariancji błędu
F = [1 dt; 0 1];                % Macierz przejścia
H = [1 0];                      % Macierz pomiaru
Q = std_u^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; % Macierz kowariancji szumu procesu dla stałej prędkości (m^2/s^4)
R = std_r^2;                    % Macierz kowariancji szumu pomiaru (m^2)

N = 300;                        % Liczba kroków czasowych

true_velocity = 100;            % (m/s)
true_pos = zeros(1, N);         % Prawdziwe położenie
turning_point = 150;            % Moment zawracania (krok czasowy)

for t = 2:N
    if t <= turning_point
        true_pos(t) = true_pos(t-1) + true_velocity * dt; 
    else
        true_pos(t) = true_pos(t-1) - true_velocity * dt; % Obiekt zawraca
    end
end

measured_pos = true_pos + std_r * randn(1, N); % szum

% Filtracja
filtered_pos = zeros(1, N);
for k = 1:N
    % Przewidywanie
    x_pred = F * x;
    P_pred = F * P * F' + Q;
    
    % Aktualizacja
    K = P_pred * H' / (H * P_pred * H' + R);
    x = x_pred + K * (measured_pos(k) - H * x_pred);
    P = (eye(size(P)) - K * H) * P_pred;
    
    % Zapis
    filtered_pos(k) = x(1);
end

% Pozycja(czas)
subplot(2, 1, 1);
plot(1:N, filtered_pos, '-r', 'LineWidth', 1.5);
hold on;; 
plot(1:N, true_pos, '--b', 'LineWidth', 1.5);
hold off;
legend('Prawdziwe Położenie', 'Zmierzone Położenie', 'Ocenione Położenie', 'Location', 'southoutside', 'Orientation', 'horizontal');
xlabel('Krok Czasowy');
ylabel('Położenie');
title('Symulacja Filtru Kalmana - estymowane położenie obiektu');

% Prędkość
re_velocity = zeros(1, N); 
re_velocity(2:end) = diff(true_pos) / dt;
estimated_velocity = diff(filtered_pos) / dt; 

% Prędkość(czas)
subplot(2, 1, 2);
plot(1:N-1, estimated_velocity, '-r', 'LineWidth', 1.5);
hold on;
plot(1:N, re_velocity, '--b', 'LineWidth', 1.5);
hold off;
legend('Oceniona Prędkość', 'Prawdziwa Prędkość', 'Location', 'southoutside', 'Orientation', 'horizontal');
xlabel('Krok Czasowy');
ylabel('Prędkość');
title('Prędkość');

sgtitle('Śledzony obiekt ze stałą prędkością, v = 100 m/s, zawracający');
