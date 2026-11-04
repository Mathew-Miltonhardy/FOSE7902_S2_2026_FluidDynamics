%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearised 2D Navier-Stokes Equation Solver for Flow Over a 
%           Stretching Sheet with Temporal Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% Show lots of decimal places
format long

% Chebyshevs and mapping parameter
N = 32; L = 2; %N was 64

% Initial Reynolds number
Re_initial = 46000;

% Reynolds Timeout
Re_timeout = 100000;

% Reynolds step-size
delR = 10; % was 10

% Start wavenumber, alpha
al_s = 0.070; % was 0.070

% End wavenumber
al_e = 0.200; % was 0.200

% Wavenumber step-size
dela = 0.0001; %was 0.0001

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Plot this stuff!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve the problem forward!
[ntcrv,e,f, ~] = StretchingSheet_LinearStability_2D_positive(al_s, al_e, dela, Re_initial, delR, Re_timeout, N,L);

% To load analysis if it has already been done
%load("./data/StabilityCurveResults")

% Prep plotting variables for first pass
figure(1); clf; hold on;
ReynoldsPlot = ntcrv(:, 2);
alphasPlot = ntcrv(:, 1);

% Remove NaNs
mask = isfinite(ReynoldsPlot);
ReynoldsPlot = ReynoldsPlot(mask);
alphasPlot = alphasPlot(mask);

% Plot
plot(ReynoldsPlot, alphasPlot, 'b', LineWidth=2)

% Make second pass to get the back half of the curve as well
% New start and end locations
newAlphaStart = max(alphasPlot);
newAlphaEnd = min(alphasPlot);
newReInitial = ReynoldsPlot(end);

% Use the reverse version that works backwards to find line
[ntcrv_reverse,e_reverse,f_reverse] = StretchingSheet_LinearStability_2D_negative(newAlphaStart, newAlphaEnd, dela, newReInitial, delR, Re_timeout, N,L);

% Prep plotting variables (ensure starts from where the other ends)
ReynoldsPlot = [newReInitial; ntcrv_reverse(:, 2)];
alphasPlot = [newAlphaStart; ntcrv_reverse(:, 1)];

% Remove NaNs
mask = isfinite(ReynoldsPlot);
ReynoldsPlot = ReynoldsPlot(mask);
alphasPlot = alphasPlot(mask);

% Set 
plot(ReynoldsPlot, alphasPlot, 'b', LineWidth=2)

% Cosmetics
xlim([0.95*Re_initial max(ReynoldsPlot)-delR])
ylim([al_s-0.05*(al_e - al_s) al_e+0.05*(al_e - al_s)])
ticksAndStuff = gca;
ticksAndStuff.FontSize = 15;
title("Neutral Stability Curve for Perturbed Flow Over a Stretching Sheet", FontSize=20, Interpreter="latex")
xlabel("Local Reynolds number, $R_\ell$", FontSize=18, Interpreter="latex")
ylabel("Streamwise wavenumber, $\alpha$", FontSize=18, Interpreter="latex")
xlim([0.95*Re_initial max(ReynoldsPlot)-delR])
ylim([al_s-0.05*(al_e - al_s) al_e+0.05*(al_e - al_s)])

% Save figure
set(figure(1), 'position', [10,50, 900,500])
saveas(figure(1), './figures/NeutralStabilityCurve', 'epsc')
saveas(figure(1), './editableFigures/NeutralStabilityCurve', 'fig')

% Save variables
save("./data/StabilityCurveResults")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Plot other stuff!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Eigenspecta %%%%%
% Find critical wave number (alpha) using min critical reynolds number
[critReynolds, indexOfCritReynolds] = min(ntcrv(:,2));
critWaveNumber = ntcrv(indexOfCritReynolds,1);

% Solve for eigenvalues
[ntcrv,e,f, ~] = StretchingSheet_LinearStability_2D_positive(critWaveNumber, critWaveNumber, dela, Re_initial, delR, Re_timeout, N,L);

% Filter infinite / NaN values
mask = isfinite(e);
e = e(mask);

% Find critical eigen value
eCritReal = ntcrv(3);
eCritIdx = find(abs(real(e)- eCritReal) < 1e-5);
eCritImag = imag(e(eCritIdx));

% Set up figure and plot
figure(2); clf; hold on;
plot(real(e),imag(e),'bx', Linewidth=1.5, HandleVisibility='off');
plot(eCritReal,eCritImag,'rx', Linewidth=1.5, DisplayName="Critical Eigenvalue");

% Plot grey line
yline(0, Color="#a3a3a3", LineStyle="--", DisplayName="Zero Line")

% Cosmetics
title("Eigenspectra Plot for Critical Wave Number", Interpreter="latex", FontSize=14)
xlabel("Real part of Eigenvalue, $\omega_r$", Interpreter="latex", FontSize=12)
ylabel("Imaginary part of Eigenvalue. $\omega_i$", Interpreter="latex", FontSize=12)
axis([0 0.16 -0.25 0.05]);
legend(Location="southwest")

% Save figure
set(figure(2), 'position', [10,50, 560,420])
saveas(figure(2), './figures/eigenspectra', 'epsc')
saveas(figure(2), './editableFigures/eigenspectra', 'fig')


%%%%% Chebyshev Mapping %%%%%
% Find Chebyshev mapping and points
for i = 1:N
    etas(i) = cos((i-1)*pi/(N-1));
    ys(i) = L * (1 - etas(i))  /(1 + etas(i));
end

% Set up figure
f = uifigure;
ax = axes(f);

% Plot main figure
plot(ax, etas, ys, 'b.', MarkerSize=10); hold on;

% Plot box where we zoomed in
rectangle(ax, 'Position',[0.5 0 0.5 25], EdgeColor='black', LineWidth=1)

% Draw dotted lines from small rectangle to new plot
rightBoxLineX = [1 0.937];
rightBoxLineY = [25 284];
leftBoxLineX = [0.5 -0.357];
leftBoxLineY = [25 284];
line(ax, rightBoxLineX,rightBoxLineY,Color="black",LineStyle="--")
line(ax, leftBoxLineX,leftBoxLineY,Color="black",LineStyle="--")

% Cosmetics
title(ax, "Chebyshev Mapping of Y-values", FontSize=14, Interpreter="latex")
xlabel(ax, "Chebyshev Parameter, $\eta$", Interpreter="latex", FontSize=12)
ylabel(ax, "Distance From Sheet, $Y$", Interpreter="latex", FontSize=12)
xlim(ax, [-1 1]) 

% Set up inner figure
inset_width = 0.5;
inset_height = 0.5;
pos = ax.Position;
inset_pos = [pos(1)+inset_width/2 pos(2)+pos(4)-inset_height/0.945 inset_width inset_height];
ax_inset = axes(f,'Position',inset_pos);

% Plot zoomed in version
plot(ax_inset, etas(1:12), ys(1:12), 'b.', MarkerSize=10);
xlim(ax_inset, [0.5 1])

% Save figure
set(f, 'position', [10,50, 560,420])
exportgraphics(f,"./figures/ChebyshevPoints.eps")
exportgraphics(f,"./editableFigures/ChebyshevPoints.fig")
saveas(f, './editableFigures/ChebyshevPoints', 'fig')


%%%%% Analytical Solutions Plots %%%%%
% Parameters
P0 = 0; % Pressure on stretching sheet

% Functions
U = @(y) exp(-y);
V = @(y) exp(-y) - 1;
P = @(y) P0 + 0.5 .* (1 - exp(-2.*y));

% Set up solutions
baseFlowys = 0:0.1:10;
Us = U(baseFlowys);
Vs = V(baseFlowys);
Ps = P(baseFlowys);

% Set up figure
figure(8); clf; plts = tiledlayout(1,3);

% Plot U solution
nexttile; hold on;
plot(Us, baseFlowys, 'b', LineWidth=2)
ticksAndStuff = gca;
ticksAndStuff.FontSize = 14;
xticks(0:0.2:1)
yticks(0:2:10)
xlabel("$x$-direction velocity, $U$", interpreter="latex", FontSize=20)
ylabel("Distance from sheet, $Y$", interpreter="latex", FontSize=20)

% Plot V solution
nexttile; hold on;
plot(Vs, baseFlowys, 'b', LineWidth=2)
ticksAndStuff = gca;
ticksAndStuff.FontSize = 14;
xticks(-1:0.2:0)
yticks(0:2:10)
xlabel("$Y$-direction velocity, $V$", interpreter="latex", FontSize=20)

% Plot P solution
nexttile; hold on;
plot(Ps, baseFlowys, 'r', LineWidth=2)
ticksAndStuff = gca;
ticksAndStuff.FontSize = 14;
xticks(0:0.1:0.5)
yticks(0:2:10)
xlabel("Pressure, $P$", interpreter="latex", FontSize=20)

% Save figure
title(plts, 'Analytical Solutions of Flow Over a Stretching Sheet', FontSize=25, interpreter="latex")
set(figure(8),'position',[20,40,1200,450])
saveas(figure(8), './figures/SplitAnalyticalSolutions', 'epsc')
saveas(figure(8), './editableFigures/SplitAnalyticalSolutions', 'fig')


%%%%% Pertubation Plots %%%%%
% Find critical wave number (alpha) using min critical reynolds number
[critReynolds, indexOfCritReynolds] = min(ntcrv(:,2));
critWaveNumber = ntcrv(indexOfCritReynolds,1);

% Solve for eigenvalues
[ntcrv,e,f, pertInd] = StretchingSheet_LinearStability_2D_positive(critWaveNumber, critWaveNumber, dela, Re_initial, delR, Re_timeout, N,L);

% Set up figure
figure(7); clf; plts = tiledlayout(1,3);
scaleValue = 1/max(abs(f(1:N,pertInd))); % scale w.r.t. largest x-momentum

% Plot scaled u pertubations
figure(7); nexttile; hold on;
uPeterb = abs(f(1:N,pertInd));
uPeterbPlot = uPeterb.*scaleValue;
plot(uPeterbPlot,ys, 'b', LineWidth=2); axis([0 1 0 10]);
ticksAndStuff = gca;
ticksAndStuff.FontSize = 14;
xticks(0:0.2:1)
yticks(0:2:10)
xlabel("$x$-direction perturbation, $u_p$", interpreter="latex", FontSize=20)
ylabel("Distance from sheet, $Y$", interpreter="latex", FontSize=20)

% Save individual plot
figure(9); clf; hold on;
plot(uPeterbPlot,ys, 'b', LineWidth=2); axis([0 1 0 10]);
ticksAndStuff = gca;
ticksAndStuff.FontSize = 14;
xticks(0:0.2:1)
yticks(0:2:10)
xlabel("$x$-direction perturbation, $u_p$", interpreter="latex", FontSize=20)
ylabel("$y$-direction (distance from plate)", interpreter="latex", FontSize=20)
set(figure(9),'position',[20,40,400,450])
exportgraphics(figure(9), './presentationFigures/uVelocityPerturbation.eps')

% Plot scaled v pertubations
figure(7); nexttile; hold on;
vPeterb = abs(f(1+N:2*N,pertInd));
vPeterbPlot = vPeterb.*scaleValue;
plot(vPeterbPlot,ys, 'b', LineWidth=2); axis([0 1 0 10]);
ticksAndStuff = gca;
ticksAndStuff.FontSize = 14;
xticks(0:0.2:1)
yticks(0:2:10)
xlabel("$Y$-direction perturbation, $v_p$", interpreter="latex", FontSize=20)

% Save individual plot
figure(10); clf; hold on;
plot(vPeterbPlot,ys, 'b', LineWidth=2); axis([0 1 0 10]);
ticksAndStuff = gca;
ticksAndStuff.FontSize = 14;
xticks(0:0.2:1)
yticks(0:2:10)
xlabel("$y$-direction perturbation, $v_p$", interpreter="latex", FontSize=20)
set(figure(10),'position',[20,40,400,450])
exportgraphics(figure(10), './presentationFigures/vVelocityPerturbation.eps')

% Plot scaled p pertubations
figure(7); nexttile; hold on;
pPeterb = abs(f(1+2*N:3*N,pertInd));
pPeterbPlot = pPeterb.*scaleValue;
plot(pPeterbPlot,ys, 'r', LineWidth=2); axis([0 1 0 10]);
ticksAndStuff = gca;
ticksAndStuff.FontSize = 14;
xticks(0:0.2:1)
yticks(0:2:10)
xlabel("Pressure perturbation, $p_p$", interpreter="latex", FontSize=20)

% Save individual plot
figure(11); clf; hold on;
plot(pPeterbPlot,ys, 'r', LineWidth=2); axis([0 1 0 10]);
xlabel("Pressure perturbation, $p_p$", interpreter="latex", FontSize=20)
ticksAndStuff = gca;
ticksAndStuff.FontSize = 14;
xticks(0:0.2:1)
yticks(0:2:10)
set(figure(11),'position',[20,40,400,450])
exportgraphics(figure(11), './presentationFigures/pPressurePerturbation.eps')

% Save figure
title(plts, 'Scaled Pertubations of Flow Over a Stretching Sheet', FontSize=25, interpreter="latex")
set(figure(7),'position',[20,40,1200,450])
saveas(figure(7), './figures/SplitAnalyticalPerturbations', 'epsc')
saveas(figure(7), './editableFigures/SplitAnalyticalPerturbations', 'fig')


%%%%% Sensitivity Analysis %%%%%
% Chebyshevs and mapping parameter
Ns = [32, 48, 64, 128];
Ls = [0.5, 1, 1.5, 2];

% % Initial Reynolds number, timeout, and step size
Re_initial_val = 48000;
Re_timeout_val = 50000;
delR_val = 10; 

% Start, end, and step size for wavenumber, alpha
al_s_val = 0.16;
al_e_val = 0.165;
dela_val = 0.0001;

% Initialise results matrix
NandLvalidaiton_wave = zeros(length(Ns), length(Ls));
NandLvalidaiton_Re = zeros(length(Ns), length(Ls));

% Solve for each combination
% For each N values
for i = 1:length(Ns)
    % Get current N
    Nc = Ns(i);

    % For each L value
    for j = 1:length(Ls)
        % Get current L
        Lc = Ls(j);

        % Solve
        [ntcrv,e,f, ~] = StretchingSheet_LinearStability_2D_positive(al_s_val, al_e_val, dela_val, Re_initial_val, delR_val, Re_timeout_val, Nc,Lc);

        % Find critical wave number (alpha) using min critical reynolds number
        [critReynolds, indexOfCritReynolds] = min(ntcrv(:,2));
        critWaveNumber = ntcrv(indexOfCritReynolds,1);

        % Write results
        NandLvalidaiton_wave(i, j) = critWaveNumber;
        NandLvalidaiton_Re(i, j) = critReynolds;
    end
end


% Function solving the problem to find the natural curve, starting at a
%    given wave number and Reynolds number, assuming stability, then
%    increasing the Reynolds number until we find instability, goes from
%    the bottom up
function [ntcrv,e,f, pertubationIndex] = StretchingSheet_LinearStability_2D_positive(al_s, al_e, dela, Re_initial, delR, Re_timeout, N,L)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Define basic variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Define complex constant
    ci = complex(0,1);
    
    % Identity matrix
    I = eye(N);
    
    % Start of Reynolds list is Re_initial 
    Res  = Re_initial;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              Discretise using Chebyshev (cheb.m)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This creates the matrix D that corresponds to the Chebyshev derivative
    [P,D] = cheb(N-1);
    % Second derivative
    D2 = D^2;
    
    % Define physical y-locations in terms of Chebyshevs u
    % Use y = L(1-eta)/(1+eta)
    for k = 1:N
        eta(k) = cos((k-1)*pi/(N-1));
        y(k) = L*(1-eta(k))./(1+eta(k));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          Base Flow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    UB = diag(exp(-y));
    DUB = diag(-exp(-y));
    VB = diag(exp(-y)-1);
    DVB = diag(-exp(-y));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    Find d eta / dy from du / dy = du / d eta * d eta / dy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DZ = -2*L./(y+L).^2;
    D2Z = 4*L./(y+L).^3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Find D = d/dy in terms of discrete Chebyshev and chain rule
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute physical Y-derivatives in terms of the 
    % Chebyshevs and transformations
    E1 = zeros(N,N);
    E2 = zeros(N,N);
    for i = 1:N
        for j = 1:N
            E1(i,j)=D(i,j)*DZ(i);
            E2(i,j)=D2(i,j)*DZ(i)*DZ(i) + D(i,j)*D2Z(i);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 Start solving this stuff!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialise progress bar
    highestPercent = 0;
    fprintf('Progress: 0%%');
    
    % Initialise neutral curve data
    ntcrv = [];
    pertubationIndex = [];
    
    % Begin for loop of wavenumber alpha
    % For each wave number from al_s to al_e by step dela
    for alpha = al_s:dela:al_e
        % Progress bar
        currPercent = round((alpha - al_s) / (al_e - al_s), 2);
        if currPercent > highestPercent
            if currPercent <= 0.1
                fprintf(repmat('\b', 1, 2));
            else
                fprintf(repmat('\b', 1, 3));
            end
            fprintf(strcat(num2str(round(100*currPercent)), '%%'))
            highestPercent = currPercent;
        end
        if alpha == al_e
            fprintf('\n')
        end
        
        % Initialise omega to just be -i 
        omega = 0-ci;  
        
        % Initialise storage of eigenvalues (omega) for each Re at current alpha
        oms = [];
        
        % Initialise Reynolds number for this loop
        Re = Res; 
    
        % Find alpha^2 for easier typing
        alpha2 = alpha^2;
    
        % flag used to update starting Reynolds for next alpha
        timeout_flag = false;
        
        % While imaginary part of eigenvalue is negative, keep checking
        % - Want imaginary part positive, unstable behaviour, so stop then
        while (imag(omega) < 0)
            % Display current Reynolds number    
            %fprintf("Reynolds Number: " + num2str(Re) + "\n")
            
            % Initialise matrices
            % - Size 3Nx3N because initially 3x3, then split into N discrete parts
            %    using Chebyshev
            A = zeros(3*N);
            B = zeros(3*N);
            
            % Set up matrices in the following form 
            % A*[u v p]' = omega*B*[u v p]'
            % Rows 1:N form the x-momentum equation 
            % Rows N+1:2N form the y-momentum equation
            % Rows 2N+1:3N form the continuity equation
            % Columns 1:N are matched to the u-velocity
            % Columns N+1:2N are matched to the v-velocity
            % Columns 2N+1:3N are matched to the pressure, p
            % E1 is D from equations, E2 is D^2 from equations
            
            %%%%%% A matrix %%%%%%
            % x-momentum equation
            A(1:N,1:N) = UB*alpha*ci + (1/Re)*(VB*E1 + UB - E2 + alpha2*I);
            A(1:N,N+1:2*N) = DUB;
            A(1:N,2*N+1:3*N) = alpha*ci*I;
        
            % A(1:N,1:N)       = ci*alpha*UB + VB*E1/Re + (alpha2*I - E2)/Re; %Check is this missing a UB?
            % A(1:N,N+1:2*N)   = DUB;
            % A(1:N,2*N+1:3*N) = ci*alpha*I;
            
            % y-momentum equation
          % A(N+1:2*N,1:N) = zeros, no need to redefine
            A(N+1:2*N,N+1:2*N) = UB*alpha*ci + (1/Re)*(VB*E1 + DVB - E2 + alpha2*I);
            A(N+1:2*N,2*N+1:3*N) = E1;
        
            % A(N+1:2*N,N+1:2*N)   = ci*alpha*UB + VB*E1/Re + DVB/Re + (alpha2*I - E2)/Re;%Again, check!
            % A(N+1:2*N,2*N+1:3*N) = E1;
            
            % Continuity equation
            A(2*N+1:3*N,1:N) = alpha*ci*I;
            A(2*N+1:3*N,N+1:2*N) = E1;
          % A(2*N+1:3*N,2*N+1:3*N) = zeros, no need to redefine
        
            % A(2*N+1:3*N,1:N)       = ci*alpha*I;
            % A(2*N+1:3*N,N+1:2*N)   = E1;
            
            %%%%%% B matrix %%%%%%
            % x-momentum equation
            B(1:N,1:N) = ci*I;
          % B(1:N,N+1:2*N) = zeros, no need to redefine
          % B(1:N,2*N+1:3*N) = zeros, no need to redefine
        
            % y-momentum equation
          % B(N+1:2*N,1:N) = zeros, no need to redefine
            B(N+1:2*N,N+1:2*N) = ci*I;
          % B(N+1:2*N,2*N+1:3*N) = zeros, no need to redefine
            
            % Continuity equation
          % B(2*N+1:3*N,1:N) = zeros, no need to redefine
          % B(2*N+1:3*N,N+1:2*N) = zeros, no need to redefine
          % B(2*N+1:3*N,2*N+1:3*N) = zeros, no need to redefine
            
        
            % Boundary conditions
            % Set rows 1,N,N+1,2N,2N+1,3N of matrices A and B to zero
            % - Rows 1 & N set u=0 on sheet and far from sheet
            % - Rows N+1 & 2N set v=0 on sheet and far from sheet
            % - Row 3N sets p=0 far from sheet
            A([1,N,N+1,2*N,2*N+1,3*N],:)=0;
            B([1,N,N+1,2*N,2*N+1,3*N],:)=0;
            
            % Include the five conditions
            % - Enforces the boundary conditions
            A(1,1) = 1; % u = 0 on the wall
            A(N,N) = 1; % u = 0 as z->infinity
            A(N+1,N+1) = 1; % v = 0 on the wall
            A(2*N,2*N) = 1; % v = 0 as z->infinity
            A(3*N,3*N) = 1; % p = 0 as z->infinity
            
            % Compute eigenvalues (e) and eigenfunction (f) of matrix equations
            [f,e] = eig(A,B);
            e = diag(e); % pull out just eigenvals, non-diagonals are 0
            
            % Calculate dominant mode - TWF wave
            % - Want to find the least negative eigenvalue
            eval = complex(0,-1.e10);
            for i = 1:3*N 
                % For each eigenvalue:
                % - If 0 < |real part| < alpha AND |imag part| < 1 (i.e. is small)
                if abs(real(e(i))) < (alpha-0.01) && abs(real(e(i))) > 0.01 && imag(e(i)) < 1
                    % - If current imag part is bigger than stored
                    if imag(e(i)) > imag(eval)
                        % Overwrite stored, save index
                        eval = e(i);
                        ind = i;
                    end
                end
            end
            
            % Largest eigenvalue (therefore will be unstable, if any of them are)
            omega = eval;
        
            % Update list of all critical eigenvalues
            oms = [oms;omega];
            
            %Update Reynolds number 
            Re = Re+delR;
            
            % If gone on too long, stop that
            if Re>Re_timeout
                timeout_flag = true;
                break
            end
        
        end %end of while loop
        
        % Determine critical conditions
        % - Check critical line exists, make NaN if not
        % - If Re times out, won't have the curve there, so NaN needed
        if length(oms) > 1
            Ren = (Re-2*delR)-imag(oms(end-1))*delR/(imag(oms(end))-imag(oms(end-1)));
            omegan = real(oms(end-1))+(Ren-(Re-2*delR))*(real(oms(end))-real(oms(end-1)))/delR;
        else
            Ren = NaN;
        end
        if timeout_flag
            Ren = NaN;
        end 
        
        % Store critical conditions for plotting
        ntcrv = [ntcrv; alpha Ren omegan];
        pertubationIndex = [pertubationIndex; ind];
        
        % Update starting Reynolds number for next alpha value
        if timeout_flag
            Res = Re_timeout - 100;
        else
            Res = max(Ren-100, Re_initial);
        end
    
    end % end of alpha for loop

    
    
end % end of function


% Function solving the problem to find the natural curve but going from the
%    top down and assuming instability and going until we find stability
function [ntcrv,e,f] = StretchingSheet_LinearStability_2D_negative(al_s, al_e, dela, Re_initial, delR, Re_timeout, N,L)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Define basic variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Define complex constant
    ci = complex(0,1);
    
    % Identity matrix
    I = eye(N);
    
    % Start of Reynolds list is Re_initial 
    Res  = Re_initial;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              Discretise using Chebyshev (cheb.m)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This creates the matrix D that corresponds to the Chebyshev derivative
    [P,D] = cheb(N-1);
    % Second derivative
    D2 = D^2;
    
    % Define physical y-locations in terms of Chebyshevs u
    % Use y = L(1-eta)/(1+eta)
    for k = 1:N
        eta(k) = cos((k-1)*pi/(N-1));
        y(k) = L*(1-eta(k))./(1+eta(k));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          Base Flow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    UB = diag(exp(-y));
    DUB = diag(-exp(-y));
    VB = diag(exp(-y)-1);
    DVB = diag(-exp(-y));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    Find d eta / dy from du / dy = du / d eta * d eta / dy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DZ = -2*L./(y+L).^2;
    D2Z = 4*L./(y+L).^3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Find D = d/dy in terms of discrete Chebyshev and chain rule
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute physical Y-derivatives in terms of the 
    % Chebyshevs and transformations
    E1 = zeros(N,N);
    E2 = zeros(N,N);
    for i = 1:N
        for j = 1:N
            E1(i,j)=D(i,j)*DZ(i);
            E2(i,j)=D2(i,j)*DZ(i)*DZ(i) + D(i,j)*D2Z(i);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 Start solving this stuff!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialise neutral curve data
    ntcrv = [];
    
    % Begin for loop of wavenumber alpha
    % For each wave number from al_s to al_e by step dela
    for alpha = al_s:-dela:al_e
        % Initialise omega to just be i 
        omega = 0+ci;  
        
        % Initialise storage of eigenvalues (omega) for each Re at current alpha
        oms = [];
        
        % Initialise Reynolds number for this loop
        Re = Res;
    
        % Find alpha^2 for easier typing
        alpha2 = alpha^2;
    
        % flag used to update starting Reynolds for next alpha
        timeout_flag = false;
        
        % While imaginary part of eigenvalue is negative, keep checking
        % - Want imaginary part negative, stable behaviour, so stop then
        while (imag(omega) > 0)
            % Display current Reynolds number    
            %fprintf("Reynolds Number: " + num2str(Re) + "\n")
            
            % Initialise matrices
            % - Size 3Nx3N because initially 3x3, then split into N discrete parts
            %    using Chebyshev
            A = zeros(3*N);
            B = zeros(3*N);
            
            % Set up matrices in the following form 
            % A*[u v p]' = omega*B*[u v p]'
            % Rows 1:N form the x-momentum equation 
            % Rows N+1:2N form the y-momentum equation
            % Rows 2N+1:3N form the continuity equation
            % Columns 1:N are matched to the u-velocity
            % Columns N+1:2N are matched to the v-velocity
            % Columns 2N+1:3N are matched to the pressure, p
            % E1 is D from equations, E2 is D^2 from equations
            
            %%%%%% A matrix %%%%%%
            % x-momentum equation
            A(1:N,1:N) = UB*alpha*ci + (1/Re)*(VB*E1 + UB - E2 + alpha2*I);
            A(1:N,N+1:2*N) = DUB;
            A(1:N,2*N+1:3*N) = alpha*ci*I;
        
            % A(1:N,1:N)       = ci*alpha*UB + VB*E1/Re + (alpha2*I - E2)/Re; %Check is this missing a UB?
            % A(1:N,N+1:2*N)   = DUB;
            % A(1:N,2*N+1:3*N) = ci*alpha*I;
            
            % y-momentum equation
          % A(N+1:2*N,1:N) = zeros, no need to redefine
            A(N+1:2*N,N+1:2*N) = UB*alpha*ci + (1/Re)*(VB*E1 + DVB - E2 + alpha2*I);
            A(N+1:2*N,2*N+1:3*N) = E1;
        
            % A(N+1:2*N,N+1:2*N)   = ci*alpha*UB + VB*E1/Re + DVB/Re + (alpha2*I - E2)/Re;%Again, check!
            % A(N+1:2*N,2*N+1:3*N) = E1;
            
            % Continuity equation
            A(2*N+1:3*N,1:N) = alpha*ci*I;
            A(2*N+1:3*N,N+1:2*N) = E1;
          % A(2*N+1:3*N,2*N+1:3*N) = zeros, no need to redefine
        
            % A(2*N+1:3*N,1:N)       = ci*alpha*I;
            % A(2*N+1:3*N,N+1:2*N)   = E1;
            
            %%%%%% B matrix %%%%%%
            % x-momentum equation
            B(1:N,1:N) = ci*I;
          % B(1:N,N+1:2*N) = zeros, no need to redefine
          % B(1:N,2*N+1:3*N) = zeros, no need to redefine
        
            % y-momentum equation
          % B(N+1:2*N,1:N) = zeros, no need to redefine
            B(N+1:2*N,N+1:2*N) = ci*I;
          % B(N+1:2*N,2*N+1:3*N) = zeros, no need to redefine
            
            % Continuity equation
          % B(2*N+1:3*N,1:N) = zeros, no need to redefine
          % B(2*N+1:3*N,N+1:2*N) = zeros, no need to redefine
          % B(2*N+1:3*N,2*N+1:3*N) = zeros, no need to redefine
            
        
            % Boundary conditions
            % Set rows 1,N,N+1,2N,2N+1,3N of matrices A and B to zero
            % - Rows 1 & N set u=0 on sheet and far from sheet
            % - Rows N+1 & 2N set v=0 on sheet and far from sheet
            % - Row 3N sets p=0 far from sheet
            A([1,N,N+1,2*N,2*N+1,3*N],:)=0;
            B([1,N,N+1,2*N,2*N+1,3*N],:)=0;
            
            % Include the five conditions
            % - Enforces the boundary conditions
            A(1,1) = 1; % u = 0 on the wall
            A(N,N) = 1; % u = 0 as z->infinity
            A(N+1,N+1) = 1; % v = 0 on the wall
            A(2*N,2*N) = 1; % v = 0 as z->infinity
            A(3*N,3*N) = 1; % p = 0 as z->infinity
            
            % Compute eigenvalues (e) and eigenfunction (f) of matrix equations
            [f,e] = eig(A,B);
            e = diag(e); % pull out just eigenvals, non-diagonals are 0
            
            % Calculate dominant mode - TWF wave
            % - Want to find the least negative eigenvalue
            eval = complex(0,-1e10);
            for i = 1:3*N 
                % For each eigenvalue:
                % - If 0 < |real part| < alpha AND |imag part| < 1 (i.e. is small)
                if abs(real(e(i))) < (alpha-0.01) && abs(real(e(i))) > 0.01 && imag(e(i)) < 1
                    % - If current imag part is bigger than stored
                    if imag(e(i)) > imag(eval)
                        % Overwrite stored, save index
                        eval = e(i);
                        ind = i;
                    end
                end
            end
            
            % Dominant eigenvalue (therefore will be stable, if any of them are)
            omega = eval;
        
            % Update list of all critical eigenvalues
            oms = [oms;omega];
            
            %Update Reynolds number 
            Re = Re+delR;
            
            % If gone on too long, stop that
            if Re>Re_timeout
                timeout_flag = true;
                break
            end
        
        end % end of while loop
        
        % Determine critical conditions
        % - Check critical line exists, make NaN if not
        % - If Re times out, won't have the curve there, so NaN needed
        if length(oms) > 1
            Ren = (Re-2*delR)-imag(oms(end-1))*delR/(imag(oms(end))-imag(oms(end-1)));
            omegan = real(oms(end-1))+(Ren-(Re-2*delR))*(real(oms(end))-real(oms(end-1)))/delR;
        else
            Ren = NaN;
        end
        if timeout_flag
            Ren = NaN;
        end 
        
        % Store critical conditions for plotting
        ntcrv = [ntcrv; alpha Ren omegan]; %did have omegan in the list at the end
        
        % Update starting Reynolds number for next alpha value
        if timeout_flag
            break % if at the end of Re window, we are done for the backwards pass
        else
            Res = max(Ren-100, Re_initial);
        end
    
    end % end of alpha for loop
    
end % end of function
