function [u1,u2,P] = RungeKutta2(K,fu1,fu2,u1,u2,N)
    global a b h
    u1(1) = u1;
    u2(1) = u2;
    % % %================
    x = a:h:b;
    % %================
%   fu1 = @(x,u1,u2) u2
%   fu2 = @(x,u1,u2) p(x)*u2 + q(x)*u1 + r(x)
    for i = 1:N
        K(1,1)  = fu1(x(i),u1(i),u2(i));
        K(1,2)  = fu2(x(i),u1(i),u2(i));
        
        K(2,1)  = fu1(x(i+1),u1(i)+h*K(1,1),u2(i)+h*K(1,2));
        K(2,2)  = fu2(x(i+1),u1(i)+h*K(1,1),u2(i)+h*K(1,2));
    
        u1(i+1) = u1(i)+(K(1,1)+K(2,1))*(h/2);
        u2(i+1) = u2(i)+(K(1,2)+K(2,2))*(h/2);
    end
    P = K;
end