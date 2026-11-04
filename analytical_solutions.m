clear all;

% Parameters
P0 = 0; % Pressure on stretching sheet

% Functions
U = @(y) exp(-y);
V = @(y) exp(-y) - 1;
P = @(y) P0 + 0.5 .* (1 - exp(-2.*y));

% Set up solutions
ys = 0:0.1:10;
Us = U(ys);
Vs = V(ys);
Ps = P(ys);

% Plot
figure(1); clf; hold on;
plot(ys, Us, 'b')
plot(ys, Vs, 'g')
plot(ys, Ps, 'r')

% Set up grid of points
[xMesh, yMesh] = meshgrid(0:0.25:3.9, 0.2:0.25:4.2);

UMesh = U(yMesh);
VMesh = V(yMesh);

% Plot flow field
figure(2); clf; hold on;
quiver(xMesh, yMesh, UMesh, VMesh, 0.5, LineWidth=1)

% Plot stretching plate
xLine = -1:0.25:4;
yLine = zeros(1, length(xLine)) + 0.05;
ULine = 0.1*ones(1, length(xLine)); %if you need arrows growing, add (xLine+0.5)*
VLine = zeros(1, length(xLine));
plot(xLine, yLine, 'k', LineWidth=8)
quiver(xLine, yLine, ULine, VLine, 0, Color='w', LineWidth=1)

% Cosmetics
ticksAndStuff = gca;
ticksAndStuff.FontSize = 15;
ylim([0 4])
xlim([0 4])
title('Analytical Solutions of Flow Over a Linear Stretching Sheet', FontSize=23, Interpreter='latex')
xlabel('$x$-direction along sheet', Interpreter='latex', FontSize=18)
ylabel('$Y$-direction perpendicular sheet', Interpreter='latex', FontSize=18)

% Savings
set(gcf,'position',[20,40,1000,450])
exportgraphics(gcf, './figures/FlowPhase.eps')


%%% Do same thing but in dimensional variables
% Parameters
xi = 1;
L = 1;
delta = 1;
% rho = 1;
% P0 = 0;

% Functions
uDim = @(x,y) xi .* L .* x .* exp(-1 / (L .* delta)) .* exp(-y);
vDim = @(x,y) xi .* L .* delta .* exp(-1 / (L .* delta)) .* exp(-y) - 1;
% pDim = @(x,y) rho .* xi.^2 .* L.^2 .* delta.^2 .* (P0 + 0.5.*(1 - exp(-1 / (L .* delta)) .* exp(-2*y)));

% Set up grid of points
[xMeshDim, yMeshDim] = meshgrid(0:0.8:10, 0.3:0.5:3);
uDimMesh = uDim(xMeshDim, yMeshDim);
vDimMesh = vDim(xMeshDim, yMeshDim);

% Pertubation things - uncomment to get random perturbations (not from
%  actual solution)
% for i = 1:13
%     for j = 1:6
%         uDimMesh(j,i) = uDimMesh(j,i) + 0.2*(j-6)*rand();
%         vDimMesh(j,i) = vDimMesh(j,i) + 0.2*(j-6)*rand();
%     end
% end

% Plot flow field
figure(3); clf; hold on;
quiver(xMeshDim, yMeshDim, uDimMesh, vDimMesh, LineWidth=2)

% Plot stretching plate
xDimLine = -1:1:10;
yDimLine = zeros(1, length(xDimLine)) + 0.035;
UDimLine = (xDimLine+0.5).*0.1.*ones(1, length(xDimLine));
VDimLine = zeros(1, length(xDimLine));
plot(xDimLine, yDimLine, 'k', LineWidth=10)
quiver(xDimLine, yDimLine, UDimLine, VDimLine, 0, Color='w', LineWidth=1.5)

% Cosmetics
ylim([0 3])
xlim([0 10])
%title('Analytical Solutions of Flow Over a Linear Stretching Sheet', FontSize=18)
xlabel('Direction along sheet', Interpreter='latex', FontSize=15)
ylabel('Direction perpendicular sheet', Interpreter='latex', FontSize=15)

% Savings
set(gcf,'position',[20,40,1000,450])
exportgraphics(gcf, './figures/DimensionalFlowPhase.eps')


% Graphs of solutions
% Parameters
P0 = 0; % Pressure on stretching sheet
xi = 1;
L = 1;
delta = 1;
rho = 1;

% Functions
U = @(y) exp(-y);
V = @(y) exp(-y) - 1;
P = @(y) P0 + 0.5 .* (1 - exp(-2.*y));

% Dimensional Version if needed
% U = @(x,y) xi .* L .* x .* exp(-1 / (L .* delta)) .* exp(-y);
% V = @(x,y) xi .* L .* delta .* exp(-1 / (L .* delta)) .* exp(-y) - 1;
% P = @(x,y) rho .* xi.^2 .* L.^2 .* delta.^2 .* (P0 + 0.5.*(1 - exp(-1 / (L .* delta)) .* exp(-2*y)));

% Set up solutions
ys = 0:0.1:6;
Us = U(ys);
Vs = V(ys);
Ps = P(ys);

% Set up figure
figure(4); clf; tiles = tiledlayout(1,3);
title(tiles, "Analytical Solutions of Flow Over a Stretching Sheet", interpreter="latex", FontSize=14)

% Plot u
nexttile; hold on; box on;
plot(Us, ys, 'b', LineWidth=2)
xlabel("$x$-direction velocity", interpreter="latex", FontSize=12)
ylabel("Distance from sheet, $Y$", interpreter="latex", FontSize=12)

% Save u
set(figure(4),'position',[20,40,400,450])
exportgraphics(figure(4), './presentationFigures/uVelocity.eps')

% Plot v
nexttile; hold on; box on;
plot(Vs, ys, 'b', LineWidth=2)
xlabel("$Y$-direction velocity", interpreter="latex", FontSize=12)

% Save v
set(figure(5),'position',[20,40,400,450])
exportgraphics(figure(5), './presentationFigures/vVelocity.eps')

% Plot p
nexttile; hold on; box on;
plot(Ps, ys, 'r', LineWidth=2)
xlabel("Pressure", interpreter="latex", FontSize=12)

% Save p
set(figure(6),'position',[20,40,400,450])
exportgraphics(figure(6), './presentationFigures/pPressure.eps')

% Save figure
set(figure(4),'position',[20,40,1200,450])
saveas(figure(4), './figures/SplitAnalyticalSolutions', 'epsc')