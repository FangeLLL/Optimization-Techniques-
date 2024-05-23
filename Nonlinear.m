clear all
close all
clc

% Define the range for x1 and x2
X1 = -5:0.01:5;
X2 = -5:0.01:5;

% Create a grid for x1 and x2
[x1, x2] = meshgrid(X1, X2);

% Define the Modified Rosenbrock function
F = 100 * (x2 - x1.^2).^2 + (6.4 * (x2 - 0.5).^2 - x1 - 0.6).^2;

% Plot the surface of the function
figure;
mesh(x1, x2, F);
title('Mesh plot of the Modified Rosenbrock function');
xlabel('x1');
ylabel('x2');
zlabel('f(x1, x2)');

% Plot the contour of the function
figure;
contourf(x1, x2, F);
title('Contour plot of the Modified Rosenbrock function');
xlabel('x1');
ylabel('x2');
hold on;


% Add initial guesses
initial_guesses = cell(3, 1);
initial_guesses{1} = -5 + 10 * rand(2, 1);  % Uniform distribution [-5, 5]
initial_guesses{2} = -5 + 10 * rand(2, 1);  % Uniform distribution [-5, 5]
initial_guesses{3} = -5 + 10 * rand(2, 1);  % Uniform distribution [-5, 5]

for i = 1:3
    plot(initial_guesses{i}(1), initial_guesses{i}(2), 'ro', 'MarkerSize', 10, 'DisplayName', ['Initial Guess ' num2str(i)]);
end

legend show;
hold off;

%% Newton-Raphson Algorithm
fprintf('Newton-Raphson Algorithm\n');
x = -5 + 10 * rand(2, 1);  % Starting point within [-5, 5]
epsilon = 10^(-4);  % Convergence tolerance

tic
fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n', x(1), x(2), func(x));
plot(x(1), x(2), 'r.');
x_next = x - inv(hessianfunc(x)) * gradfunc(x);
fprintf('k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n', x_next(1), x_next(2), func(x_next), norm(gradfunc(x_next)));
plot(x_next(1), x_next(2), 'r*');
k = 3;

while (norm(gradfunc(x_next)) > epsilon)
    x = x_next;
    x_next = x - inv(hessianfunc(x)) * gradfunc(x);
    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, error=%f\n', k, x_next(1), x_next(2), func(x_next), norm(gradfunc(x_next)));
    plot(x_next(1), x_next(2), 'r*');
    k = k + 1;
end

toc
title('Newton-Raphson Algorithm');
set(gca, 'fontsize', 35);
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Hestenes-Stiefel Algorithm
figure;
contourf(x1, x2, F);
hold on;

fprintf('Hestenes-Stiefel Algorithm\n');
x = -5+10*rand(2, 1);  % Starting point
epsilon = 10^(-4);  % Convergence tolerance

tic
fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n', x(1), x(2), func(x));
plot(x(1), x(2), 'r.');
g = gradfunc(x);
d = -g;

% alpha argmin procedure
alpha_vals = 0:0.01:1;
funcalpha = zeros(length(alpha_vals), 1);
for i = 1:length(alpha_vals)
    funcalpha(i) = func(x + alpha_vals(i) * d);
end
[val, ind] = min(funcalpha);
alpha = alpha_vals(ind);
% end of alpha argmin procedure

x_next = x + alpha * d;
g_next = gradfunc(x_next);
beta = (g_next' * (g_next - g)) / (d' * (g_next - g));
d_next = -g_next + beta * d;

fprintf('k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n', x_next(1), x_next(2), func(x_next), norm(gradfunc(x_next)));
plot(x_next(1), x_next(2), 'r*');
k = 3;

while (norm(gradfunc(x_next)) > epsilon)
    x = x_next;
    g = g_next;
    d = d_next;

    % alpha argmin procedure
    alpha_vals = 0:0.01:1;
    funcalpha = zeros(length(alpha_vals), 1);
    for i = 1:length(alpha_vals)
        funcalpha(i) = func(x + alpha_vals(i) * d);
    end
    [val, ind] = min(funcalpha);
    alpha = alpha_vals(ind);
    % end of alpha argmin procedure

    x_next = x + alpha * d;
    g_next = gradfunc(x_next);
    beta = (g_next' * (g_next - g)) / (d' * (g_next - g));
    d_next = -g_next + beta * d;

    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, error=%f\n', k, x_next(1), x_next(2), func(x_next), norm(gradfunc(x_next)));
    plot(x_next(1), x_next(2), 'r*');
    k = k + 1;
end

toc
title('Hestenes-Stiefel Algorithm');
set(gca, 'fontsize', 35);
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Polak-Ribiere Algorithm
figure;
contourf(x1, x2, F);
hold on;

fprintf('Polak-Ribiere Algorithm\n');
x = -5 + 10 * rand(2, 1);  % Starting point within [-5, 5]
epsilon = 10^(-4);  % Convergence tolerance

tic
fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n', x(1), x(2), func(x));
plot(x(1), x(2), 'r.');
g = gradfunc(x);
d = -g;

% Armijo line search parameters
sigma = 0.1;
beta = 0.7;
alpha = 1;

k = 2;

while (true)
    % Armijo line search
    alpha = 1;  % Start with alpha = 1
    while func(x + alpha * d) > func(x) + sigma * alpha * g' * d
        alpha = beta * alpha;  % Reduce alpha
    end

    x_next = x + alpha * d;
    g_next = gradfunc(x_next);

    % Compute the Polak-Ribiere beta
    beta_pr = (g_next' * (g_next - g)) / (g' * g);
    
    % Update the direction
    d = -g_next + beta_pr * d;

    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n', k, x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)));
    plot(x_next(1), x_next(2), 'r*');

    % Check for termination
    if norm(g_next) <= epsilon || abs(func(x_next) - func(x)) <= epsilon
        break;
    end

    % Update variables for next iteration
    x = x_next;
    g = g_next;
    k = k + 1;
end

toc
title('Polak-Ribiere Algorithm');
set(gca, 'fontsize', 35);
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Fletcher-Reeves Algorithm
figure;
contourf(x1, x2, F);
hold on;

fprintf('Fletcher-Reeves\n');
x = -5 + 10 * rand(2, 1);  % Starting point within [-5, 5]
epsilon = 10^(-4);  % Convergence tolerance

tic
fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n', x(1), x(2), func(x));
plot(x(1), x(2), 'r.');
g = gradfunc(x);
d = -g;

% Armijo line search parameters
sigma = 0.1;
beta = 0.7;
alpha = 1;

k = 2;

while (true)
    % Armijo line search
    alpha = 1;  % Start with alpha = 1
    while func(x + alpha * d) > func(x) + sigma * alpha * g' * d
        alpha = beta * alpha;  % Reduce alpha
    end

    x_next = x + alpha * d;
    g_next = gradfunc(x_next);

    % Compute the Fletcher-Reeves beta
    beta_fr = (g_next' * g_next) / (g' * g);
    
    % Update the direction
    d = -g_next + beta_fr * d;

    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n', k, x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)));
    plot(x_next(1), x_next(2), 'r*');

    % Check for termination
    if norm(g_next) <= epsilon || abs(func(x_next) - func(x)) <= epsilon
        break;
    end

    % Update variables for next iteration
    x = x_next;
    g = g_next;
    k = k + 1;
end

toc
title('Fletcher-Reeves Algorithm');
set(gca, 'fontsize', 35);
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

