clc
clear
close
epsilon = 0.1;
x = -5:0.01:5;
N = length(x);
for i = 1:N
if x(i)<-1
    y(i) = -1;
    if x(i)>=-1-epsilon
        z(i) = -(-x(i) + (-x(i)-1)/epsilon - ((-x(i))^2-1)/(2*epsilon));
    elseif x(i)< -1-epsilon
        z(i) = -1-epsilon/2;
    end
elseif (x(i)>=-1) && (x(i)<=1)
    y(i) = x(i);
    z(i) = x(i);
elseif x(i)>1
    y(i) = 1;
    if x(i)<=1+epsilon
        z(i) = x(i) + (x(i)-1)/epsilon - (x(i)^2-1)/(2*epsilon);
    elseif x(i)>=1+epsilon
        z(i) = 1 + epsilon/2;
    end
end
end
plot(x,y,'r-',x,z,'b-.','linewidth',1.5); grid on;
legend('saturation','modified')