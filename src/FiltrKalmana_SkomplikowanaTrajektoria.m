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

true_velocity = zeros(1, N);    % (m/s)
change_point1 = 100;            % punkt zmiany 1
change_point2 = 200;            % punkt zmiany 2

% prędkość zmieniająca się w czasie
for i = 1:N
    if i <= change_point1
        true_velocity(i) = 25;
    elseif i <= change_point2
        true_velocity(i) = 100;
    else
        true_velocity(i) = 0.1;
    end
end

true_pos = zeros(1, N);         % Prawdziwe położenie
turning_point = 150;            % Moment zawracania (krok czasowy)

% położenie
for i = 2:N
    if i <= turning_point
        true_pos(i) = true_pos(i-1) + true_velocity(i-1) * dt;
    else
        true_pos(i) = true_pos(i-1) - true_velocity(i-1) * dt;
    end
end

measured_pos = true_pos + std_r * randn(1, N); % szum


% Filtracja
filtered_pos = zeros(1, N);
filtered_velocity = zeros(1, N);
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

% Prędkość
re_velocity = zeros(1, N); 
re_velocity(2:end) = diff(true_pos) / dt; 
estimated_velocity = zeros(1, N); 
estimated_velocity(2:end) = diff(filtered_pos) / dt; 

% Wykres pozycji
subplot(2, 1, 1);
plot(1:N, filtered_pos, '-r', 'LineWidth', 1.5);
hold on;
plot(1:N, true_pos, '--b', 'LineWidth', 1.5);
hold off;
legend('Estymowane Położenie', 'Prawdziwe Położenie', 'Location', 'southoutside', 'Orientation', 'horizontal');
xlabel('Krok Czasowy');
ylabel('Położenie');
title('Symulacja Filtru Kalmana - estymowane położenie obiektu');

% Prędkość(czas)
subplot(2, 1, 2);
plot(1:N, estimated_velocity, '-r', 'LineWidth', 1.5);
hold on;
plot(1:N, re_velocity, '--b', 'LineWidth', 1.5);
hold off;
legend('Estymowana Prędkość', 'Prawdziwa Prędkość', 'Location', 'southoutside', 'Orientation', 'horizontal');
xlabel('Krok Czasowy');
ylabel('Prędkość');
title('Symulacja Filtru Kalmana - estymowana prędkość obiektu');

sgtitle('Śledzony obiekt z zmienną prędkością i zawracający w połowie');
