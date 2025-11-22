L = 32;   % <<<--- SET YOUR L HERE

% your nonlinear function BEFORE rescaling
g = @(x) -0.5*(x./(L-1)).^2 + x./(L-1);

% compute min/max of g over 0..L
xs = linspace(0, L, 2000);
gmin = min(g(xs));
gmax = max(g(xs));

% rescaled function from 0 to 1
f = @(x) (g(x) - gmin) ./ (gmax - gmin);

% % ========== INVERSE FUNCTION ==========
% % For y in [0,1], find x such that f(x) = y using fzero
% 
% finv = @(y) arrayfun(@(yy) fzero(@(xx) f(xx) - yy, ...
%                  L * yy), ...  % initial guess scales with y
%                  y);
% 
% % ========== PLOT THE INVERSE ==========
% yvals = linspace(0,1,400);
% xvals = finv(yvals);
% 
% plot(yvals, xvals, 'LineWidth', 1.5)
% 
% hold on
% plot from x = 0 to x = L, y in [0,1]
fplot(f, [0 L])
ylim([0 1])
grid on
