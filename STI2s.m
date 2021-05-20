function z=STI2s(n,lambda,sigma,zeta,c,mu,beta,um,us,un,T, init,p,P,y)
%STI2s(n,lambda,sigma,f,c,mu,beta,um,us,un,T, init,p,P,y) simulates the evolution of a continuous-time behavioural SIS epidemic model with a 2-stage negotiation process
%n is the number of nodes
% lambda is the per-contact infection probability
% sigma is the complaisance
% zeta is the risk perception
% c is the cost
% mu is the detection rate
% beta is the recovery rate
% um is the control effort in marketing
% us is the control effort in screening
% un is the control effort in partner notification
% T is the duration of the simulation
% init is the initial fraction of infected individuals
% p is the faithfulness
% P are the partners
% y is a variable such that if set to 'y'the function produces the plot
% z is the long-run epidemic prevalence

c=c-um;
if length(init)==1
    x=zeros(1,n);
    x(randperm(n,round(n*init)))=2;
else
    x=init;
end
if length(sigma)==1
    sigma=sigma*ones(n,1);
end
if length(p)==1
    temp=p;
for i=1:n
if P(i)~=i
p(i)=temp;
else 
    p(i)=0;
end
end
end
zie(1)=sum(x)/2;
zsp(1)=0;
zit(1)=0;
zse(1)=n-zie(1);
time(1)=0;
t=1;
energy=(mu+us+c+beta+zeta+1+un)*n;
R=[zeta us c mu beta un];
R=sommacumulativa(R)/sum(R);
dp=n/energy;
while time(t)<T
    t=t+1;
    zse(t)=zse(t-1);
    zsp(t)=zsp(t-1);
    zie(t)=zie(t-1);
    zit(t)=zit(t-1);
    time(t)=time(t-1)+randexp(energy);
    if rand<=dp %contact  %SE=0, SP=1, IE=2, IP=3
        i=randi(n);
if rand<=p(i) %partner
            j=P(i);
            flag=0;
        else %random
         j=modn(i+randi(n-1),n); %select a different node
         flag=1;
        end
        if rand>p(j) || flag==0
            state_i=x(i);
        state_j=x(j);
        if state_i==0 && state_j==1 
            if flag*rand<=sigma(j)
            x(j)=0;
            zsp(t)=zsp(t)-1;
            zse(t)=zse(t)+1;
            elseif flag*rand<=sigma(i) %counterproposal
            x(i)=1;
            zsp(t)=zsp(t)+1;
            zse(t)=zse(t)-1;
            end
        end
         if state_i==0 && state_j==2 && rand<=lambda
            x(i)=2;
            zie(t)=zie(t)+1;
            zse(t)=zse(t)-1;
         end
          if state_i==0 && state_j==3 
              if flag*rand<=sigma(j)
                if rand<=lambda
                    x(i)=2;
                    x(j)=2;
                    zie(t)=zie(t)+2;
                    zse(t)=zse(t)-1;
                else
                    x(j)=2;
                    zie(t)=zie(t)+1;
                end
              elseif flag*rand<=sigma(i) %counterproposal
                  x(i)=1;
                  zsp(t)=zsp(t)+1;
                  zse(t)=zse(t)-1;
              end
          end
           if state_i==1 && state_j==0 
               if flag*rand<=sigma(j)
                x(j)=1;
                zsp(t)=zsp(t)+1;
                zse(t)=zse(t)-1;
               elseif flag*rand<=sigma(i)
               x(i)=0;
                zsp(t)=zsp(t)-1;
                zse(t)=zse(t)+1;
               end
           end
            if state_i==1 && state_j==2 
                if flag*rand<=sigma(j)
            x(j)=3;
            zie(t)=zie(t)-1;
                elseif flag*rand<=sigma(i)
                    x(i)=2;
                    zie(t)=zie(t)+1;
                    zsp(t)=zsp(t)-1;
                end
            end
            if state_i==2 && state_j==0 && rand<=lambda
            x(j)=2;
            zie(t)=zie(t)+1;
            zse(t)=zse(t)-1;
            end
            if state_i==2 && state_j==3 
                if flag*rand<=sigma(j)
            x(j)=2;
            zie(t)=zie(t)+1;
                elseif flag*rand<=sigma(i)
                    x(i)=3;
                    zie(t)=zie(t)-1;
                end
            end
            if state_i==2 && state_j==1 
                if flag*rand<=sigma(j)
              if rand<=lambda
             x(j)=2;
             zie(t)=zie(t)+1;
             zsp(t)=zsp(t)-1;
              else
             x(j)=0;
             zse(t)=zse(t)+1;
             zsp(t)=zsp(t)-1;
              end
                elseif flag*rand<=sigma(i)
                    x(i)=3;
                    zie(t)=zie(t)-1; 
                end
            end
          if state_i==3 && state_j==0 
              if flag*rand<=sigma(j)
            x(j)=1;
            zse(t)=zse(t)-1;
            zsp(t)=zsp(t)+1;
              elseif flag*rand<=sigma(i)
                  if rand<=lambda
                      x(i)=2;
                      x(j)=2;
                      zse(t)=zse(t)-1;
                      zie(t)=zie(t)+2;
                  else
                      x(i)=2;
                      zie(t)=zie(t)+1;
                  end
              end
          end
            if state_i==3 && state_j==2 
                if flag*rand<=sigma(j)
            x(j)=3;
            zie(t)=zie(t)-1;
                elseif flag*rand<=sigma(i)
                    x(i)=2;
                    zie(t)=zie(t)+1;
                end
            end
        end
    else %recover
        r=randvett(R);
        switch r
            case 1 %zeta
                i=randi(n);
                if x(i)==0 && rand<=zit(t)/n
                    x(i)=1;
                    zse(t)=zse(t-1)-1;
                    zsp(t)=zsp(t-1)+1;
                elseif x(i)==2 && rand<=zit(t)/n
                    x(i)=3;
                    zie(t)=zie(t-1)-1;
                end
            case 2 %us
                i=randi(n);
                switch x(i)
                    case 0
                        x(i)=1;
                        zse(t)=zse(t-1)-1;
                        zsp(t)=zsp(t-1)+1;
                    case 2
                        x(i)=4;
                        zie(t)=zie(t-1)-1;
                        zit(t)=zit(t-1)+1;
                    case 3
                        x(i)=4;
                        zit(t)=zit(t-1)+1;
                end
            case 3 %c 
                i=randi(n);
                if x(i)==1
                    x(i)=0;
                    zse(t)=zse(t-1)+1;
                    zsp(t)=zsp(t-1)-1;
                elseif x(i)==3
                    x(i)=2;
                    zie(t)=zie(t-1)+1;
                end
             case 4 %mu
                i=randi(n);
                switch x(i)
                    case 2
                        x(i)=4;
                        zie(t)=zie(t-1)-1;
                        zit(t)=zit(t-1)+1;
                    case 3
                        x(i)=4;
                        zit(t)=zit(t-1)+1;
                end
            case 5 %beta
                 i=randi(n);
                if x(i)==4
                    x(i)=1;
                    zit(t)=zit(t-1)-1;
                    zsp(t)=zsp(t-1)+1;
                end 
            otherwise %un
                i=randi(n);
                if x(P(i))==4
                 switch x(i)
                    case 0
                        x(i)=1;
                        zse(t)=zse(t-1)-1;
                        zsp(t)=zsp(t-1)+1;
                    case 2
                        x(i)=4;
                        zie(t)=zie(t-1)-1;
                        zit(t)=zit(t-1)+1;
                    case 3
                        x(i)=4;
                        zit(t)=zit(t-1)+1;
                 end
                end
        end  
    end
end
 [zse,t]=reducev2(zse/n,time,400);
  [zsp,t]=reducev2(zsp/n,time,400);
   [zie,t]=reducev2(zie/n,time,400);
      [zit,t]=reducev2(zit/n,time,400);
      if y=='y'
  figure
  plot(t,zse,'b')
  hold on
  plot(t,zsp,'g')
    plot(t,zie,'r')
      plot(t,zit,'m')
    plot(t,1-zse-zsp-zie-zit,'y')
      end
      z=1-zse(end)-zsp(end);
end

