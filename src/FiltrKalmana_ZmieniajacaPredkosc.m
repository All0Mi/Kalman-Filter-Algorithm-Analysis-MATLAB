std_r = 10;                     % odchylenie pomiaru (z dokładnością do 10m)
std_u = 1;                      % odchylenie procesu (m/s)
dt = 1;                         % krok czasowy (s)
x = [0; 0];                     % Początkowa estymata stanu [pozycja; prędkość] (m, m/s)
P = eye(2);                     % Początkowa estymata kowariancji błędu
F = eye(2);                     % Macierz przejścia
H = [1 0];                      % Macierz pomiaru
Q = std_u^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; % Macierz kowariancji szumu procesu dla stałej prędkości (m^2/s^4)
R = std_r^2;                    % Macierz kowariancji szumu pomiaru (m^2)

N = 300;                        % Liczba kroków czasowych


true_pos = zeros(1, N);          
true_velocity = zeros(1, N);    

% zmieniająca się prędkość
t = linspace(0, 2*pi, N);
true_velocity = 20 + 100*sin(t)+2*cos(3*t);   

for i = 2:N
    true_pos(i) = true_pos(i-1) + true_velocity(i-1) * dt; 
end

% Szum
measured_pos = true_pos + std_r * randn(1, N); % szum

% Inicjalizacja
filtered_pos = zeros(1, N);
x_pred = x;
P_pred = P;

% Filtracja
for k = 1:N
    % Aktualizacja 
    F = [1 dt; 0 1];    
    F(1, 2) = true_velocity(k);     
    
    % Przewidywanie
    x_pred = F * x;
    P_pred = F * P * F' + Q;
    
    % Aktualizacja
    K = P_pred * H' / (H * P_pred * H' + R);
    x = x_pred + K * (measured_pos(k) - H * x_pred);
    P = (eye(size(P)) - K * H) * P_pred;
    
    filtered_pos(k) = x(1);
end

% Pozycja(czas)
subplot(2, 1, 1);
plot(1:N, filtered_pos, '-r', 'LineWidth', 1.5);
hold on;
plot(1:N, true_pos, '--b', 'LineWidth', 1.5);
hold off;
legend('Zmierzone Położenie', 'Prawdziwe Położenie', 'Location', 'southoutside', 'Orientation', 'horizontal');
xlabel('Krok Czasowy');
ylabel('Położenie');
title('Symulacja Filtru Kalmana - estymowane położenie obiektu');

% Prędkość
estimated_velocity = diff(filtered_pos) / dt;

% Prędkość(czas)
subplot(2, 1, 2);
plot(1:N-1, estimated_velocity, '-r', 'LineWidth', 1.5);
hold on;
plot(1:N, true_velocity, '--b', 'LineWidth', 1.5);
hold off;
legend('Oceniona Prędkość', 'Prawdziwa Prędkość', 'Location', 'southoutside', 'Orientation', 'horizontal');
xlabel('Krok Czasowy');
ylabel('Prędkość');
title('Prędkość');

sgtitle('Śledzony obiekt ze zmieniającą się prędkością');
