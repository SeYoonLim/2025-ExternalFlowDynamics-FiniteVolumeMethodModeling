clc;
clear;
close all;

% Input values
a = 3;
b = 5;
h = 4;

% Compute trapezoid area directly
A = (a + b) * h / 2;

% Print result
fprintf('a = %.2f\n', a);
fprintf('b = %.2f\n', b);
fprintf('h = %.2f\n', h);
fprintf('Trapezoid area A = %.2f\n', A);

%% Visualization 1: sine wave
x = linspace(0, 2*pi, 300);
y = sin(x);

figure;
plot(x, y, 'LineWidth', 1.8);
grid on;
xlabel('x');
ylabel('sin(x)');
title('Sine Wave Example');

%% Visualization 2: trapezoid shape
x_trap = [0, b, (b + a)/2, (b - a)/2, 0];
y_trap = [0, 0, h, h, 0];

figure;
fill(x_trap, y_trap, [0.7 0.85 1.0], 'EdgeColor', 'k', 'LineWidth', 1.5);
grid on;
axis equal;
xlabel('x');
ylabel('y');
title(sprintf('Trapezoid Shape (Area = %.2f)', A));

text(b/2, h/2, sprintf('A = %.2f', A), ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');