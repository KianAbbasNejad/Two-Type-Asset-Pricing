% Heterogeneous Expectations Asset Pricing Model 
% This code calculates the time series, attractor and bifurcation diagram
% of a two type asset pricing model as described in Hommes (2013, Ch.6) 
% and his lectures at University of Amsterdam in winter of 2019.
% (Written by Kian Abbas Nejad)

clearvars
clc
close all

tmax=10000;
T = 1:tmax; % time 
X = zeros(1,tmax); % Deviation from RE Equilibrium
U1 = zeros(1,tmax); %Utility of fundamentalists
U2 = zeros(1,tmax); % utility of trend followers
N1 = zeros(1,tmax); %fraction of fundamentalists
N2 = zeros(1,tmax); % fraction of trend followers

% Parameters
Beta = 3.6; 
R = 1.1;
g1 = 0;
b1 = 0; 
g2 = 1.2;
b2 = 0; 
C1 = 1;
C2= 0;
mu = 0;
sigma = 0.01;

%Initial Conditions -- need 3 boundary conds
X(1:2)=[0.1 0.2];
N1(1:3)=0.2;

% Formulae
fu1= @(t,x) (x(t)-R*x(t-1))*(g1*x(t-2)+b1-R*x(t-1))-C1;
fu2= @(t,x) (x(t)-R*x(t-1))*(g2*x(t-2)+b2-R*x(t-1))-C2;
fn1 = @(t,u1,u2,beta) exp(beta*u1(t-1))/(exp(beta*u1(t-1))+exp(beta*u2(t-1)));
fn2 = @(t,n1) 1-n1(t);
fx = @(t,n1,n2,x) (n1(t)*(g1*x(t-1)+b1) + n2(t)*(g2*x(t-1)+b2) + normrnd(mu,sigma))/R;

% ODE System Starts at period 3
for i=3:tmax
    N2(i) = fn2(i,N1);
    X(i)=fx(i,N1,N2,X);
    U1(i)= fu1(i,X);
    U2(i)= fu2(i,X);
    N1(i+1)= fn1(i+1,U1,U2,Beta);
end
N1(end)=[]; % Remove the last unwanted iteration


% TIME SERIES 
fig1 = figure(1);
fig1.Name='Time Series';
plot(X(1:400),'k')
grid on
title('Time Series','FontSize',14,'interpreter','latex');
ylabel('$x_t$','FontSize',16,'interpreter','latex');
xlabel('$t$','FontSize',16,'interpreter','latex');

fig2 = figure(2);
scatter(X,N1,1,'k')
fig2.Name='Attractor';
grid on
title('Attractor (With Small Noise)','FontSize',14,'interpreter','latex');
ylabel('$n_{1,t}$','FontSize',16,'interpreter','latex');
xlabel('$x_t$','FontSize',16,'interpreter','latex');


% BIFURCATION: 
Beta = 2:0.002:4;


% initialise bifurcation diagram
fig3=figure(3);
fig2.Name = 'Bifurcation';
title('Bifurcation Diagram','FontSize',14,'interpreter','latex');
ylabel('$x_t$','FontSize',14,'interpreter','latex');
xlabel('$\beta$','FontSize',14,'interpreter','latex');

hold on
for k=1:length(Beta) % first loop iterates through parameter
    for i=3:1000 % same as before
        N2(i) = fn2(i,N1);
        X(i)=fx(i,N1,N2,X);
        U1(i)= fu1(i,X);
        U2(i)= fu2(i,X);
        N1(i+1)= fn1(i+1,U1,U2,Beta(k));
    end
    N1(end)=[]; % Remove the last unwanted iteration
    % create a vector of Beta(k) in order to plot it:
    B = ones(1,501)*Beta(k);
    scatter(B,X(500:1000),0.03,'k');
end



