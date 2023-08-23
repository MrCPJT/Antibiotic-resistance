function [ym,zm] = monoplot(t0,dt,t1,mu1,mu2,M11,M22)

te = t0:dt:t1;

ym = zeros(t1/dt,1); ym(1) = (1);  % S. Aureus population
zm = zeros(t1/dt,1); zm(1) = (1);  % P. Aeruginosa population

for i = 1:(numel(te)-1)

    ym(i+1) = ym(i) + dt * ym(i)*mu1 * ( 1 - ym(i)/(M11) );
    zm(i+1) = zm(i) + dt * zm(i)*mu2 * ( 1 - zm(i)/(M22) );
    
end

end