% Damped simple pendulum

% Equation of motion: 
% d^2theta/dt^2 + 2*gamma*omega0*dtheta/dt + omega0^2*sin(theta)=0

% Rewrite as two, first order ODEs:
% dy1/dt = y2
% dy2/dt= -2*gamma*omega0*y2-omega0_sq*sin(y1)

function [period,sol] = damped_oscillator(omega0, gamma, theta0, thetadot0)

m = 1;
g = 9.81;
R = g/omega0^2;
T_small_angle = 2*pi/omega0;
b = 2*m*gamma*omega0;

y0 = [theta0, thetadot0];                   % Initial condition
N = 10;                                     % Number of cycles
steps = 500;
tspan = linspace(0,N*T_small_angle,steps);

opts = odeset('refine',6);
[t,y]=ode45(@f,tspan,y0,opts,gamma,omega0);        % Solve ODE
sol = [t,y];

ind1 = y(:,2).*circshift(y(:,2), [-1 0]) <= 0;         % Require thetadot=0
ind2 = abs(y(:,1)) >= abs(0.01*y(1,1));                % Require amplitude at least 1% of theta0
lgcl_ind = ind1 & ind2;                                % Logical indicator
period = 2*mean(diff(t(lgcl_ind)));                    % Find t satisifying indicator

%Plot theta-time
figure
plot(t,y(:,1),'b',t,y(:,2),'r','linewidth',2);
legend('\theta','d\theta/dt')
title(['\theta and d\theta/dt v.s. time with \gamma = ' num2str(gamma)])
ylabel('\theta or d\theta/dt')
xlabel('t(s)')

function dydt=f(t,y,gamma,omega0)
dydt = [y(2);-2*gamma*omega0*y(2)-omega0^2*sin(y(1))];
