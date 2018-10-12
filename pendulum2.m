function [period,sol] = pendulum2(omega,theta0,thetad0,grph) 
% Finds the period of a nonlinear pendulum given the length of the pendulum
% arm and initial conditions. All angles in radians.

%Setting initial conditions
if nargin==0
    error('Must input length and initial conditions')
end
if nargin==1
   theta0 = pi/2;
   thetad0=0;
   grph=0;
end
if nargin==2
    thetad0 = 0;
    grph=1;
end
if nargin==3
    grph=1;
end
if nargin==4
    grph=0;
end
m = 1;
g = 9.81;
R = g/omega^2;
%omega = sqrt(g/R);

T = 2*pi/omega;
% number of oscillations to graph
N = 10;


tspan = [0 N*T];
%opts = odeset('events',@events,'refine',6); %Here for future event finder
opts = odeset('refine',6);
r0 = [theta0 thetad0];
[t,w] = ode45(@proj,tspan,r0,opts,g,R);
sol = [t,w];

lgcl_ind = w(:,2).*circshift(w(:,2), [-1 0]) <= 0;    % Logical indicator
period= 2*mean(diff(t(lgcl_ind)));                    % Find t satisifying indicator

if nargin==3
    % Calculate total energy
    Ek = 1/2 * m * (R*w(:,2)).^2;
    Ep = m*g * R * (1-cos(w(:,1)));
    E = Ek + Ep;
end

if grph % Plot all the graphs
    
    % Plot total energy v.s. time
    figure
    plot(t,E)
    title('Total energy v.s. time')
    xlabel('t')
    ylabel('Total energy')

    
end
end
%-------------------------------------------
%
function rdot = proj(t,r,g,R)
    rdot = [r(2); -g/R*sin(r(1))];
end