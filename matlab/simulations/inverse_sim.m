%--------------------------------------
% Visualization of f(x), f^{-1}(x), and f^{-1}(f(x)) = x
%--------------------------------------

% Define the domain
x = linspace(0,1,400);

% Example function: choose something smooth and invertible on [0,1]
% You may modify f if you want a different shape.
% f = @(x) x.^3 + 0.1*sin(4*pi*x);   % strictly increasing on [0,1]
f = @(x) x.^3 + 0.1*x.^2;
% Numerical inverse using interpolation
y = f(x);
f_inv = @(yq) interp1(y, x, yq, 'linear');

% Compute inverse values
x2 = f_inv(f(x));

% Plot
figure; hold on; grid on; box on;

% Plot f(x)
plot(x, f(x), 'LineWidth', 2, 'Color', [0.1 0.4 0.9]);

% Plot f^{-1}(x)
plot(x, f_inv(x), 'LineWidth', 2, 'Color', [0.9 0.4 0.1]);

% Plot identity line
plot(x, x, '--k', 'LineWidth', 1.5);

xlabel('Input');
ylabel('Output');
title('Visualization of f(x), f^{-1}(x), and the Identity f^{-1}(f(x)) = x');
legend('f(x)', 'f^{-1}(x)', 'Identity: x', 'Location', 'NorthWest');

axis([0 1 0 1]);  % restrict view to [0,1] × [0,1]
