%==========================================================================
% Tutorial Modeling and Simulation
% Topic#1: Quantiles simulation
% Authors: M.Trifonov 
% Email: trifonov.m@yahoo.com
% Date(dd-mm-yyyy): 25-03-2023
%==========================================================================
clc, clear all, close all
% initial coordinates
xl = -4.1;
xr = 4.2; %0.9908

% parameters of the normal distribution
h=0.01;
i=0;
mu = 0;
sigma = 1;
X  = -10:0.01:10;
Yp = normpdf(X, mu, sigma);
% quantile
fi=0.1285;

% equations of the line
x1 = 0.5:0.01:10;
x2 = -10:0.01:0.5;
y1 = lin_func1(x1);
y2 = lin_func2(x2);
% boundaries coordinates
f1 = lin_func2(xl);
f2 = lin_func1(xr);

% Case 1
if f1 > fi & f2 > fi

    figure(1);hold on;
    plot([-10 10],[fi fi],'c')
    plot(X,Yp),grid on; 
    plot(x1,y1,'b',x2,y2,'r');grid on;
    xlabel('random value');ylabel('probability density');

    while fi <= f1 & fi <= f2  
        f1 = lin_func2(xl);
        f2 = lin_func1(xr);
        i = i + 1;

        if f1 > f2 % Case 1.1
            xl = xl + h;
            f1 = lin_func2(xl);
        else  % Case 2.2    
            xr = xr - h;
            f2 = lin_func1(xr);
        end
        f1i(i) = f1;
        f2i(i) = f2;
        counter(i) = i;

        P = normcdf([xl xr]);  
        P_res(i) = P(2)-P(1);

        drawnow;
        figure(1);plot(xl,f1,'b*',xr,f2,'r*'),hold on,grid on
    end
    % convergence of the probability
    figure(2)
    plot(counter, P_res,'b');hold on;grid on; 
    xlabel('Number of steps'); ylabel('Probability');
    % convergence of the estimate fi
    figure(3);hold on;
    plot(counter, f1i, 'b');grid on;
    plot(counter, f2i, 'r');grid on;

    P_res  
end

%Case 2
if f1<=fi & f2<=fi
    
    figure(1);hold on;
    plot([-10 10],[fi fi],'c')
    plot(X,Yp),grid on; 
    plot(x1,y1,'b',x2,y2,'r');grid on;
    xlabel('random value');ylabel('probability density');
    
    while fi >= f1 & fi >= f2  
        f1 = lin_func2(xl);
        f2 = lin_func1(xr);
        i = i + 1;
    
        if f1 < f2 % Case 2.1
            xl = xl - h;
            f1 =  lin_func2(xl);
        else % Case 2.2
            xr = xr + h;
            f2 = lin_func1(xr);
        end
        f1i(i) = f1;
        f2i(i) = f2;
        counter(i) = i;

        P = normcdf([xl xr]);  
        P_res(i) = P(2)-P(1);
        
        drawnow;
        figure(1);plot(xl,f1,'b*',xr,f2,'r*'),hold on,grid on

    end
    % convergence of the probability
    figure(2)
    plot(counter, P_res,'b');hold on;grid on; 
    xlabel('Колличество шагов'); ylabel('Вероятность');
    % convergence of the estimate fi
    figure(3);hold on;
    plot(counter, f1i, 'b');grid on;
    plot(counter, f2i, 'r');grid on;

    P_res  
end

function [ y ] = lin_func1(x)
k1 = 0.05;
b  = -0.025;
y = k1*x + b;
end

function [ y ] = lin_func2(x)
k1 = -0.048;
b  = 0.02399;
y = k1*x + b;
end
