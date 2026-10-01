clear;
clc;

%Tarea 1: Cargar y graficar
data = readtable('datos_motor.txt');

%Extraer datos
t = data{:,2};
u = data{:,3};
y = data{:,4};

%Graficar entrada y salida
figure;
plot(t, u, 'b', t, y, 'r', 'LineWidth', 1.5);
xlabel('Tiempo [s]');
ylabel('Amplitud');
title('Respuesta del sistema');
legend('u(t)','y(t)', 'Location', 'best');
grid on;
hold on;

%Tarea 2: Lineas y tangente

%Detectar el escalon
idx_step = find(u > u(1), 1);
t_step = t(idx_step);

%Promediar antes del escalon
y0 = mean(y(1:idx_step-1));

%Promediar las ultimas 20 muestras
y_inf = mean(y(end-19:end));

%Graficar las lineas
yline(y0, '--', 'Linea base', ...
    'Color', [1.0 0.8 0.2]);

yline(y_inf, '--', 'Linea del 100%', ...
    'Color', [1.0 0.8 0.2]);

%Seleccionar puntos de la tangente
disp('Ingresa 2 puntos en la seccion de mas pendiente');
[tp, yp] = ginput(2);

%Calcular pendiente e intercepto
m = (yp(2) - yp(1)) / (tp(2) - tp(1));
b = yp(1) - m * tp(1);

%Graficar la tangente
y_tang = m * t + b;
plot(t, y_tang, 'm', 'LineWidth', 1.2);

%Calcular los cortes
t_corte_base = (y0 - b) / m;
t_corte_100 = (y_inf - b) / m;

%Marcar los cortes
plot(t_corte_base, y0, 'ko', 'MarkerFaceColor', 'k');
plot(t_corte_100, y_inf, 'ko', 'MarkerFaceColor', 'k');

%Mostrar datos de la Tabla 1.4
disp('--- Valores para la Tabla 1.4 ---');
fprintf('t1 = %.4f s, y1 = %.4f\n', tp(1), yp(1));
fprintf('t2 = %.4f s, y2 = %.4f\n', tp(2), yp(2));
fprintf('Pendiente (m) = %.4f\n', m);
fprintf('Corte con linea base (t_A o theta) = %.4f s\n', t_corte_base);
fprintf('Corte con linea 100%% (t_B) = %.4f s\n', t_corte_100);

%Ajustar los ejes
ylim([-0.2, 1.2]);
xlim([0, max(t)]);

%Tarea 3: Identificar y comparar

%Calcular la ganancia
delta_y = y_inf - y0;
delta_u = u(end) - u(1);
K = delta_y / delta_u;

%Calcular los niveles
y_28 = y0 + 0.283 * delta_y;
y_63 = y0 + 0.632 * delta_y;

%Encontrar los cruces
idx_28 = find(y >= y_28, 1);
t_28 = t(idx_28);

idx_63 = find(y >= y_63, 1);
t_63 = t(idx_63);




%t_A y t_B: cortes de la tangente en tiempo absoluto
%t1 = t_28 - t_step; t2 = t_63 - t_step
%Metodo Ziegler & Nichols
%theta = t_A - t_step
theta_zn = t_corte_base - t_step;
%tau = t_B - t_A
tau_zn = t_corte_100 - t_corte_base;
%G(s) = K*exp(-theta*s)/(tau*s + 1)
G_zn = tf(K, [tau_zn 1], 'InputDelay', theta_zn);
%Y(s) = G(s)*U(s)
ym_zn = lsim(G_zn, u, t);
%ECM = (1/N)*sum((y_k - y_hat_k)^2)
ecm_zn = mean((y - ym_zn).^2);

%Metodo Miller
%theta = t_A - t_step
theta_miller = t_corte_base - t_step;
%tau = t_63.2 - theta; t_63.2 medido desde el escalon
tau_miller = (t_63 - t_step) - theta_miller;
%G(s) = K*exp(-theta*s)/(tau*s + 1)
G_miller = tf(K, [tau_miller 1], 'InputDelay', theta_miller);
%Y(s) = G(s)*U(s)
ym_miller = lsim(G_miller, u, t);
%ECM = (1/N)*sum((y_k - y_hat_k)^2)
ecm_miller = mean((y - ym_miller).^2);

%Metodo analitico
%tau = 1.5*(t2 - t1)
tau_analitico = 1.5 * (t_63 - t_28);
%theta = t2 - tau
theta_analitico = (t_63 - t_step) - tau_analitico;
%G(s) = K*exp(-theta*s)/(tau*s + 1)
G_analitico = tf(K, [tau_analitico 1], 'InputDelay', theta_analitico);
%Y(s) = G(s)*U(s)
ym_analitico = lsim(G_analitico, u, t);
%ECM = (1/N)*sum((y_k - y_hat_k)^2)
ecm_analitico = mean((y - ym_analitico).^2);



%Graficar los modelos y los datos
figure;
plot(t, u, 'b', t, y, 'r', 'LineWidth', 1.5);
hold on;
plot(t, ym_zn, 'y--', 'LineWidth', 1.2);
plot(t, ym_miller, 'm-.', 'LineWidth', 1.2);
plot(t, ym_analitico, 'g:', 'LineWidth', 1.5);

yline(y0, 'w--', 'Linea base');
yline(y_inf, 'w--', 'Linea del 100%');

ylim([-0.1, max(y)*1.2]);
xlim([0, max(t)]);
grid on;

title('Comparacion de modelos FOTD vs Datos Reales');
xlabel('Tiempo [s]');
ylabel('Amplitud');
legend('Escalon u(t)', 'Respuesta real y(t)', ...
    'Ziegler-Nichols', 'Miller', 'Analitico', 'Location', 'best');

%Mostrar parametros
disp('--- Valores para la Tabla 1.5 (Parametros) ---');
fprintf('K (comun) = %.4f\n', K);
fprintf('Ziegler & Nichols: theta = %.4f s, tau = %.4f s\n', theta_zn, tau_zn);
fprintf('Miller:            theta = %.4f s, tau = %.4f s\n', theta_miller, tau_miller);
fprintf('Analitico:         theta = %.4f s, tau = %.4f s\n', theta_analitico, tau_analitico);

%Mostrar errores
disp('--- Valores para la Tabla 1.6 (ECM) ---');
fprintf('ECM Ziegler & Nichols = %.6f\n', ecm_zn);
fprintf('ECM Miller            = %.6f\n', ecm_miller);
fprintf('ECM Analitico         = %.6f\n', ecm_analitico);

%Datos para Simulink
datos_reales = timeseries(y, t);