clear;
close all;

% Parametri
A=0.2; % Ampiezza per u(t)
L=3;
n=3;

% Funzione u(s)
u=@(s) A*sin(pi*n/2/L*s);

% Integrandi
integrand_x=@(s) sin(u(s));
integrand_y=@(s) cos(u(s));

% Intervallo di t
t_values=linspace(0,L,200); % Genera 200 punti da 0 a L

x_values = zeros(size(t_values)); % Pre-allocazione per efficienza
y_values = zeros(size(t_values)); % Pre-allocazione per efficienza

% Calcolo degli integrali per ogni valore di t
for i=1:length(t_values)
    t_val=t_values(i);
    
    % Integra da 0 a t_val
    x_t=integral(integrand_x,0,t_val);
    y_t=integral(integrand_y,0,t_val);
    
    x_values(i)=x_t;
    y_values(i)=y_t;
end

% Salva i dati in un file (es. curve_data.dat)
% Ogni riga conterrà una coppia (x, y) separata da spazio
fileID = fopen('data3.dat', 'w');
for i = 1:length(x_values)
    fprintf(fileID, '%.6f %.6f\n', x_values(i), y_values(i));
end
fclose(fileID);

disp('Dati della curva salvati in curve_data.dat');

% (Opzionale) Visualizza la curva in MATLAB
figure;
plot(x_values, y_values, 'b-', 'LineWidth', 1.5);
pbaspect([1,1,1]);
daspect([1,1,1]);
xlabel('x');
ylabel('y');
title('Curva Parametrizzata');