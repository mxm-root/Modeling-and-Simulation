%==========================================================================
% Tutorial Modeling and Simulation
% Example #3: Statistical analysis of space ship trajectories
% Authors: M.Trifonov
% Email: trifonov.m@yahoo.com
% Date(dd-mm-yyyy): 09-06-2023
%==========================================================================
clc;clear;close all;

syms x1 x2 x3 x4  f1 f2 f3 f4  g  f  x

% jacobian
f1=-g*sin(x2);
f2=-g/x1*cos(x2);
f3=x1*sin(x2);
f4=x1*cos(x2);
f=[f1;f2;f3;f4];
x=[x1;x2;x3;x4];
A=jacobian(f,x)

clear all

% initial mx, Kx: random initial state vector
mx=zeros(4,1); Kx=zeros(4); 
mx(1,1)=500; % m/s
mx(2,1)=-0.2; % rad
mx(3,1)=5000; % m
mx(4,1)=0; % m
Kx(1,1)= 100; % m/s 
Kx(2,2)=0.01; % rad^2
Kx(3,3)=10000; % m^2
Kx(4,4)=0; %  m^2
g=9.81; % m/s^2
dt=1; % step time
tk=5000; % end time
j=0; % counter

for t = 0:dt:tk % time loop
    
    j = j + 1;
    
    f(1,1)=-g*sin(mx(2,1));
    f(2,1)=-g/mx(1,1)*cos(mx(2,1));
    f(3,1)=mx(1,1)*sin(mx(2,1));
    f(4,1)=mx(1,1)*cos(mx(2,1));
    
    % state matrix A
    A=[ 0, -g*cos(mx(2,1)), 0, 0;...
         g/mx(1,1)^2*cos(mx(2,1)), g/mx(1,1)*sin(mx(2,1)), 0, 0;...
         sin(mx(2,1)), mx(1,1)*cos(mx(2,1)), 0, 0;...
         cos(mx(2,1)), -mx(1,1)*sin(mx(2,1)), 0, 0]; % матрица А
     
    mx = mx+f*dt; % the Euler integration method for Kx, mx
    
    % check altitude (must be > 0)
    if mx(3,1) < 0
        warning('Stop condition: altitude average less than 0');
        break
    end
    
    time(j) = t; 
    
    mV(j) = mx(1,1); % plotting mean(V)
    mL(j) = mx(4,1); % plotting mean(L)
    
    Kx = Kx+(A*Kx+Kx*A')*dt; % the Euler integration method for Kx
    sigV(j) = sqrt(Kx(1,1)); % plotting sigma(V)
    sigL(j) = sqrt(Kx(4,4)); % plotting sigma(L)
    
end
% plotting
figure(1);
plot(time,mV,time,mV+3*sigV,time,mV-3*sigV);grid on;
legend('\itm_V', '\itm_V \pm 3\sigma_V');
xlabel('Time (s)');ylabel('\it Airspeed V (m)')
figure(2);
plot(time,mL,time,mL+3*sigL,time,mL-3*sigL);grid on;hold on;
legend('\itm_L', '\itm_L \pm 3\sigma_L');
xlabel('Time (s)');ylabel('\it Range L (m)')
