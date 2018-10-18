% Forced damped simple pendulum

% Equation of motion: 
% d^2theta/dt^2 + 2*gamma*dtheta/dt + omega0^2*sin(theta)= A0 cos(wt)

% Rewrite as two, first order ODEs:
% dy1/dt = y2
% dy2/dt= -2*gamma*y2-omega0_sq*sin(y1) + A0*cos(wt)

function [period,sol,A_steady,sol_steady] = forced_oscillator(...
                        omega0, gamma, A0, w, theta0, thetadot0, grph)

if nargin<=4
    error('Must input initial conditions')
end
if nargin==5
    thetadot0=0;
    grph=1;
end
if nargin==6
    grph=1;
end
if nargin==7
    grph=0;
end

m = 1;
g = 9.81;
R = g/omega0^2;
T_small_angle = 2*pi/omega0;
b = 2*m*gamma;

y0 = [theta0, thetadot0];                   % Initial condition
N = 40;                                     % Number of cycles
steps = 500;
tspan = linspace(0,N*T_small_angle,steps);

if gamma<=3
    opts = odeset('refine',6);
    [t,y]=ode45(@f,tspan,y0,opts,omega0,gamma,A0,w);        % Solve ODE
else
    opts = odeset('events',@events,'refine',6);
    [t,y]=ode45(@f,tspan,y0,opts,omega0,gamma,A0,w);        % Stop integration when equilibrium
end
sol = [t,y];

if w==0
    ind = y(:,2).*circshift(y(:,2), [-1 0]) <= 0;     % Require thetadot=0
    period = 2*mean(diff(t(ind)));                    % Find t satisifying indicator
    A_steady = A0;
    sol_steady = sol;
else
    t_steady = t(250:end);
    y_steady = y(250:end,:);
    sol_steady = [t_steady,y_steady];
    ind = y_steady(:,2).*circshift(y_steady(:,2), [-1 0]) <= 0;     % Require thetadot=0
    period = 2*mean(diff(t_steady(ind)));
    A = abs(y_steady(ind,1));
    A = A(1:end-2);
    A_steady = mean(A);
end

if grph
    figure
    plot(t,y(:,1),'b',t,y(:,2),'r','linewidth',2);
    title(['\theta and d\theta/dt v.s. t with \gamma = ' num2str(gamma) ', A0 = 1, \omega = ' num2str(w)])
    ylabel('\theta')
    xlabel('t(s)')
end

function dydt=f(t,y,omega0,gamma,A0,w)
dydt = [y(2);-2*gamma*y(2)-omega0^2*sin(y(1))+A0*cos(w*t)];


function [value,isterminal,dir] = events(t,y,omega0,gamma,A0,w)
% Locate the time when angle passes through 0.0001 in a 
% decreasing direction and stop integration.
value = y(1)-0.0001;    % Detect theta = 0.0001
isterminal = 1;         % Stop the integration
dir = -1;               % Negative direction only
