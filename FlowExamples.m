%%%%% Plotting rotating flow %%%%%
% Parameter value
L = 0.5;

% Define velocity
xp = @(x,y) -L*y;
yp = @(x,y) L*x;

% Set up grid of points
[xMesh, yMesh] = meshgrid(-6:1:6, -6:1:6);
xpMesh = xp(xMesh, yMesh);
ypMesh = yp(xMesh, yMesh);

% Plot
figure(4); clf; hold on;
quiver(xMesh, yMesh, xpMesh, ypMesh, LineWidth=2)

% Cosmetics
ticksAndStuff = gca;
ticksAndStuff.FontSize = 15;
xticks(-5:5)
yticks(-5:5)
title('Analytical Solutions of a Rotating Flow', FontSize=25, Interpreter='latex')
xlabel('$x$-direction', Interpreter='latex', FontSize=20)
ylabel('$y$-direction', Interpreter='latex', FontSize=20)
plot(0, 0, 'b.')
ylim([-5.1 5.1])
xlim([-5.1 5.1])

% Save figure
set(figure(4), 'position', [10,50, 900,500])
saveas(figure(4), './figures/RotatingFlowPhase', 'epsc')


%%%%% Plotting stagnation flow %%%%%
% Parameter value
a = 2;

% Define velocity
xp = @(x,y) a*x;
yp = @(x,y) -a*y;

% Set up grid of points
[xMesh, yMesh] = meshgrid(-5:1:5, -5:1:5);
xpMesh = xp(xMesh, yMesh);
ypMesh = yp(xMesh, yMesh);

% Plot
figure(5); clf; hold on;
quiver(xMesh, yMesh, xpMesh, ypMesh, LineWidth=2)

% Cosmetics
ticksAndStuff = gca;
ticksAndStuff.FontSize = 15;
xticks(-5:5)
yticks(-5:5)
title('Analytical Solutions of a Stagnant Point Flow', FontSize=25, Interpreter='latex')
xlabel('$x$-direction', Interpreter='latex', FontSize=20)
ylabel('$y$-direction', Interpreter='latex', FontSize=20)
plot(0, 0, 'b.')
ylim([-5 5])
xlim([-5 5])

% Save figure
set(figure(5), 'position', [10,50, 900,500])
saveas(figure(5), './figures/StagnantFlowPhase', 'epsc')
