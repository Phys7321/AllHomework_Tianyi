function [period,sol] = pendulum1(omega,theta0,thetad0,grph) 
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

T= 2*pi/omega;
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
    Delta_n = (E(:) - E(1)) / E(1);

    idx = find(lgcl_ind);
    idx = idx(1:2:end);     % Find the starting index of every cycle
    % Average energy in every cycle
    EK = [];
    EP = [];
    for i=1:N
        EK(i) = mean(Ek(idx(i):idx(i+1)));
        EP(i) = mean(Ep(idx(i):idx(i+1)));
    end
end

if grph % Plot all the graphs
    
    % Plot total energy v.s. time
    figure
    scatter(t,Delta_n)
    title('Relative change in total energy in first oscillation')
    xlim([0,T])
    xlabel('t')
    ylabel('\Delta')
    
    % Plot position and velocity v.s. time
    figure
    subplot(2,1,1)
    plot(t,w(:,1),'b')
    title('Position v.s. time')
    xlabel('t')
    ylabel('\theta')
    subplot(2,1,2)
    plot(t,w(:,2),'r')
    title('Velocity v.s. time')
    xlabel('t')
    ylabel('d\theta / dt')
    
    % Plot average KE and KP
    figure
    plot(1:N,EK,'ro',1:N,EP,'bo',1:N,EK+EP,'ko-')
    title('Average EK and EP during complete cycle')
    legend('EK_{avg}','EP_{avg}','Etotal_{avg}','Location','best')
    xlabel('Number of cycle')
    ylabel('Energy')
    
end
end
%-------------------------------------------
%
function rdot = proj(t,r,g,R)
    rdot = [r(2); -g/R*sin(r(1))];
end