function P = partners_net(n,p)
%  n number of nodes, p, probability of being in a relationship
P=1:n;
np=round(n*p/2);
temp=randperm(n,2*np);
for i=1:np
P(temp(1+2*(i-1)))=temp(2*i);
P(temp(2*i))=temp(1+2*(i-1));
end
end