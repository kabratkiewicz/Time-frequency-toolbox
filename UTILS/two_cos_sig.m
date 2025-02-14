function [a, if1,if2,s1,s2,s] = two_cos_sig(t)

N = length(t);
a = 1.7-gausswin(N,1.6);
a = a(:);
t = t(:);

if1 = 125*t + 3*sin(2*pi.*t.* 3 - pi/4); 
if2 = -125*t + 2*cos(2*pi.*t * 8 + pi/4); 

s1 = a.*exp(2*pi*1i.*(if1));
s2 = a.*exp(2*pi*1i.*(if2));

s = s1+s2;
end