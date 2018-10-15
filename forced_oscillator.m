% Forced damped simple pendulum

% Equation of motion: 
% d^2theta/dt^2 + 2*gamma*dtheta/dt + omega0^2*sin(theta)= A0 cos(wt)

% Rewrite as two, first order ODEs:
% dy1/dt = y2
% dy2/dt= -2*gamma*y2-omega0_sq*sin(y1) + A0*cos(wt)

function [period,sol] = forced_oscillator(omega0, gamma, A0, w, theta0, thetadot0, grph)

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
N = 15;                                     % Number of cycles
steps = 500;
tspan = linspace(0,N*T_small_angle,steps);

if gamma<=3
    opts = odeset('refine',6);
    [t,y]=ode45(@f,tspan,y0,opts,omega0,gamma,A0,w);        % Solve ODE
else
    opts = odeset('events',@events,'refine',6);
    [t,y]=ode45(@f,tspan,y0,opts,omega0,gamma,A0,w);        % Solve ODE
end

sol = [t,y];

%if w==0
    ind1 = y(:,2).*circshift(y(:,2), [-1 0]) <= 0;         % Require thetadot=0
    ind2 = abs(y(:,1)) >= 0.0001;                          % Require amplitude at least 0.0001
    lgcl_ind = ind1 & ind2;                                % Logical indicator
    period = 2*mean(diff(t(lgcl_ind)));                    % Find t satisifying indicator
%else
%    opts = odeset('refine',6);
%    [t_unpert,y_unpert]=ode45(@f,tspan,y0,opts,omega0,gamma,0,0);  % Find the steady state through                                                   
%    ind3 = abs(y_unpert(:,1)) <= 0.0001;                           %   unperturbed    
%    ind1 = y(:,2).*circshift(y(:,2), [-1 0]) <= 0;                 %   case
%    lgcl_ind = ind1 & ind3;                                    
%    period = 2*mean(diff(t(lgcl_ind)));  
%end

if grph
    figure
    plot(t,y(:,1),'b','linewidth',2);
    title(['\theta v.s. t with \gamma = 0.5, A0 = 1, \omega = ' num2str(w)])
    ylabel('\theta')
    xlabel('t(s)')
end

function dydt=f(t,y,omega0,gamma,A0,w)
dydt = [y(2);-2*gamma*y(2)-omega0^2*sin(y(1))+A0*cos(w*t)];


function [value,isterminal,dir] = events(t,y,omega0,gamma,A0,w)
% Locate the time when height passes through zero in a 
% decreasing direction and stop integration.
value = y(1)-0.0001;   % Detect height = 0
isterminal = 1;   % Stop the integration
dir = -1;   % Negative direction only
