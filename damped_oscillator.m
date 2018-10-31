% Damped simple pendulum

% Equation of motion: 
% d^2theta/dt^2 + gamma*dtheta/dt + omega0^2*theta=0

% Rewrite as two, first order ODEs:
% dy1/dt = y2
% dy2/dt= -gamma*y2-omega0_sq*y1

function [period,sol] = damped_oscillator(omega0, gamma, theta0, thetadot0, grph)

if nargin==0||nargin==1
    error('Must input length and initial conditions')
end
if nargin==2
   theta0 = 1;
   thetadot0=0;
   grph=0;
end
if nargin==3
    thetadot0 = 0;
    grph=1;
end
if nargin==4
    grph=1;
end
if nargin==5
    grph=0;
end

m = 1;
g = 9.81;
R = g/omega0^2;
T_small_angle = 2*pi/omega0;

y0 = [theta0, thetadot0];                   % Initial condition
N = 10;                                     % Number of cycles
steps = 500;
tspan = linspace(0,N*T_small_angle,steps);


%------Solve ODE and find period-------
if gamma<6
    opts = odeset('refine',6);
    [t,y]=ode45(@f,tspan,y0,opts,gamma,omega0);        % Solve ODE
else
    opts = odeset('events',@events,'refine',6);
    [t,y]=ode45(@f,tspan,y0,opts,gamma,omega0);        % Solve ODE
end

sol = [t,y];

ind = y(:,2).*circshift(y(:,2), [-1 0]) <= 0;         % Require thetadot=0
period = 2*mean(diff(t(ind)));                        % Find t satisifying indicator


%------Calculate Energy-------
Ek = 1/2 * m * (R*y(:,2)).^2;
Ep = m*g * R * (1-cos(y(:,1)));
E = Ek + Ep;
N_cycles = round((N*T_small_angle)/period);
Start=1;
E_avg = zeros(N_cycles+1,1);
t_avg = zeros(N_cycles+1,1);
E_avg(1) = E(1);
t_avg(1) = t(1);
for i=1:N_cycles
    remainder=mod(steps,N_cycles);
    step=round((steps-remainder)/N_cycles);
    E_avg(i+1) = mean(E(Start:(Start+step-1)));
    t_avg(i+1) = mean(t(Start:(Start+step-1)));
    Start=Start+step;
end


if grph
    figure
    plot(t,y(:,1),'b',t,y(:,2),'r','linewidth',2);
    legend('\theta','d\theta/dt')
    title(['\theta and d\theta/dt v.s. time with \gamma = ' num2str(gamma)])
    ylabel('\theta or d\theta/dt')
    xlabel('t(s)')
    
    figure
    plot(t,E,'b',t_avg,E_avg,'ro--');
    title(['Total Energy v.s. time with \gamma = ' num2str(gamma)])
    ylabel('Energy')
    xlabel('t(s)')
end

function dydt=f(t,y,gamma,omega0)
dydt = [y(2);-gamma*y(2)-omega0^2*y(1)];


function [value,isterminal,dir] = events(t,y,gamma,omega0)
% Locate the time when angle passes through 0.0001 in a 
% decreasing direction and stop integration.
value = y(1)-0.0001;    % Detect theta = 0.0001
isterminal = 1;         % Stop the integration
dir = -1;               % Negative direction only
