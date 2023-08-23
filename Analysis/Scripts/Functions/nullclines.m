function [x,y,u,v,intersection,w1] = nullclines(mu1,mu2,M11,M22,M12,M21)

% Function to visualise the nullclines of our gLV system
% 

x1 = linspace(-1000000,1000000,10000001);

w1 = M12*(1-x1/M11);
w2 = M22*(1-x1/M21);

inter=find(abs(w1-w2)<min(abs(w1-w2)+1));

try 
    intersection=inter(1);
catch
    intersection=1;
end

f = @(t,Y) [mu1*Y(1)*(1-Y(1)/M11-Y(2)/M12);
mu2*Y(2)*(1-Y(2)/M22-Y(1)/M21)];

y1 = linspace(x1(intersection)-50,x1(intersection)+50,15);
y2 = linspace(w1(intersection)-50,w1(intersection)+50,15);

[x,y] = meshgrid(y1,y2);
size(x)
size(y)

u = zeros(size(x));
v = zeros(size(x));

t=0; 

for i = 1:numel(x)

    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);

end


end