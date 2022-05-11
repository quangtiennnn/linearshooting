function [u1,u2,P] = Euler(fu1,fu2,u1,u2,N)
    global a b h
    u1(1) = u1;
    u2(1) = u2;
    % % %================
    x = a:h:b;
    % %================
%   fu1 = @(x,u1,u2) u2
%   fu2 = @(x,u1,u2) p(x)*u2 + q(x)*u1 + r(x)
    for i = 1:N
        u1(i+1)  = u1(i) + h*fu1(x(i),u1(i),u2(i));
        u2(i+1)  = u2(i) + h*fu2(x(i),u1(i),u2(i));
    end
end