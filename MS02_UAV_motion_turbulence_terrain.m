%==========================================================================
% Tutorial Modeling and Simulation
% Example#2: An UAV motion taking into account turbulence and terrain
% Authors: M.Trifonov 
% Email: trifonov.m@yahoo.com
% Date(dd-mm-yyyy): 16-03-2019
%==========================================================================
clear; clc;
% parameters of the PID controller
Kp=3.159; Kd=4.4552; Ki=-0.0024; 
% parameters of the UAV dynamic model
K=0.2;T=0.15;
E=0.5;V=280;
A=40;  
% parameters of the terrain
Krel=300; Trel=10; Erel=5;
% parameters of the turbulence
Ot=12; L=400; l=L/V;
% parameter of the sensor
Nn=0.1;
% initial data
h0 = 60;htr = 25;
% simulation
N = 100; % number of simulations

for i = 1:N
    
    a1 = 23341 + i;
    a2 = 24341 + i;
    a3 = 25341 + i;
    % start of the simulink model
    sim('MS02_SimModel_01.slx');
    % save of the simulation data
    H(:,i) = hout(:,2);
    
end

% calculation of mean and standard deviation
for j = 1:length(H(:,1))
    
    mx(j) = mean(H(j,:));
    dx(j) = cov(H(j,:));
    
end

% plotting
tout = hout(:,1);
x_vektor = [tout',fliplr(tout')];
y_vektor = [mx+3*sqrt(dx),fliplr(mx-3*sqrt(dx))];
figure(1); hold on; grid on;
patch = fill(x_vektor,y_vektor, [243 169 114]./255);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.4);
plot(tout,mx,'m', 'LineWidth', 2);
plot(tout,hout(:,2),'b:', 'LineWidth', 0.5);
xlabel('Time [s]'); ylabel('Flight altitude [m]')