clc
clear all

global a b N h;
N = 10;
a = 1;
b = 2;
h = (b-a)/N

syms x u1 u2 u3 u4;
K = zeros(4,2);
p(x) = -2./x;
q(x) = 2./(x.^2);
r(x) = sin(log(x))./(x.^2);
alpha = 1;
beta = 2;
    %R1<=x<=R2
    %y(2) = V1;
    %y(4) = 0;
% %y Exact
c2 =  (1/70)*(8-12*sin(log(2))-4*cos(log(2)))
c1 = 11/10 - c2
yex = @(x) c1.*x+ c2./(x.^2) - (3./10).*sin(log(x))-(1./10).*cos(log(x));


%4 - IVP: 
fu1 = @(x,u1,u2) u2
fu2 = @(x,u1,u2) p(x).*u2 + q(x).*u1 + r(x)
    u1(1) = alpha;
    u2(1) = 0;
fu3 = @(x,u3,u4) u4 
fu4 = @(x,u3,u4) p(x).*u4 + q(x).*u3
    u3(1) = 0;
    u4(1) = 1;
% % %================
x = a:h:b;
% % %================
%Runge-Kutta4
[u1 u2 K] = RungeKutta4(K,fu1,fu2,u1(1),u2(1));
[u3 u4 K] = RungeKutta4(K,fu3,fu4,u3(1),u4(1));
%Runge-Kutta2
% [u1 u2 K] = RungeKutta2(K,fu1,fu2,u1(1),u2(1));
% [u3 u4 K] = RungeKutta2(K,fu3,fu4,u3(1),u4(1));
%Euler
% [u1 u2] = Euler(fu1,fu2,u1(1),u2(1));
% [u3 u4] = Euler(fu3,fu4,u3(1),u4(1));

double(u1) 
double(u2);
double(u3)
double(u4);
% % %================
w1t = alpha;
w2t = (beta - u1(N+1))/u3(N+1);
% % %================
w1 = u1 + w2t*u3;
w2 = u2 + w2t*u4;
double(w1);

y = yex(x);
for j = 2:N
    E(j-1) = abs((w1(j)-y(j)));
end
double(E)