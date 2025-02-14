function [a, if1,if2,s1,s2,s] = two_chirps_sig(t)

N = length(t);
a = 0.6+gausswin(N,1.6);
a = a(:);
t = t(:);

if1 = 20*t + 200.*t.^2./2; 
if2 = -175*t + 100.*t.^2./2; 

s1 = a.*exp(2*pi*1i.*(if1));
s2 = a.*exp(2*pi*1i.*(if2));

s = s1+s2;
end
