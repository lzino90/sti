function [z,time]=STI_adnc(n,ar,lambda,sigma,f,ua,c,mu,beta,T, init)
%SIS_adnc(n,a,eta,mu,lambda, init) follows the evolution of a continuous-time SIS epidemics
%on a discrete-distrubution ADN with n nodes divided in activation classes
%with rates a. eta is the fraction of nodes in each class. lambda the
%per-contact infection probability and mu the recovery rate. init is the
%initial condition, it can be the fraction x_1(0) or the initial configuration x(0). T is the total observation time.
if length(init)==1
    x=zeros(1,n);
    x(randperm(n,round(n*init)))=2;
else
    x=init;
end
if length(ar)==1
    ar=ar*ones(1,n);
end
zie(1)=sum(x)/2;
zsp(1)=0;
zit(1)=0;
zse(1)=n-zie(1);
time(1)=0;
t=1;
energy=(mu+ua+c+beta+f)*n+sum(ar);
R=[f ua c mu beta];
R=sommacumulativa(R)/sum(R);
dp=sum(ar)/energy;
sc=sommacumulativa(ar)/sum(ar);
while time(t)<T
    t=t+1;
    zse(t)=zse(t-1);
    zsp(t)=zsp(t-1);
    zie(t)=zie(t-1);
    zit(t)=zit(t-1);
    time(t)=time(t-1)+randexp(energy);
    if rand<dp %contact  %SE=0, SP=1, IE=2, IP=3
        i=randvett(sc);
        j=randi(n);
        state_i=x(i);
        state_j=x(j);
        if state_i==0 && state_j==1 && rand<sigma
            x(j)=0;
            zsp(t)=zsp(t)-1;
            zse(t)=zse(t)+1;
        end
         if state_i==0 && state_j==2 && rand<lambda
            x(i)=2;
            zie(t)=zie(t)+1;
            zse(t)=zse(t)-1;
         end
          if state_i==0 && state_j==3 && rand<sigma
              if rand<lambda
             x(i)=2;
             x(j)=2;
             zie(t)=zie(t)+2;
             zse(t)=zse(t)-1;
              else
             x(j)=2;
             zie(t)=zie(t)+1;
              end
          end
           if state_i==1 && state_j==0 && rand<sigma
            x(j)=1;
            zsp(t)=zsp(t)+1;
            zse(t)=zse(t)-1;
           end
            if state_i==1 && state_j==2 && rand<sigma
            x(j)=3;
            zie(t)=zie(t)-1;
            end
            if state_i==2 && state_j==0 && rand<lambda
            x(j)=2;
            zie(t)=zie(t)+1;
            zse(t)=zse(t)-1;
            end
            if state_i==2 && state_j==3 && rand<sigma
            x(j)=2;
            zie(t)=zie(t)+1;
            end
            if state_i==2 && state_j==1 && rand<sigma
              if rand<lambda
             x(j)=2;
             zie(t)=zie(t)+1;
             zsp(t)=zsp(t)-1;
              else
             x(j)=0;
             zse(t)=zse(t)+1;
             zsp(t)=zsp(t)-1;
              end
            end
          if state_i==3 && state_j==0 && rand<sigma
            x(j)=1;
            zse(t)=zse(t)-1;
            zsp(t)=zsp(t)+1;
          end
            if state_i==3 && state_j==2 && rand<sigma
            x(j)=3;
            zie(t)=zie(t)-1;
            end
    else %recover
        r=randvett(R);
        switch r
            case 1 %f 
                i=randi(n);
                if x(i)==0 && rand<zit(t)/n
                    x(i)=1;
                    zse(t)=zse(t-1)-1;
                    zsp(t)=zsp(t-1)+1;
                elseif x(i)==2 && rand<zit(t)/n
                    x(i)=3;
                    zie(t)=zie(t-1)-1;
                end
            case 2 %ua
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
            case 3 %c (to be done)
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
            otherwise
                i=randi(n);
                if x(i)==4
                    x(i)=1;
                    zit(t)=zit(t-1)-1;
                    zsp(t)=zsp(t-1)+1;
                end
        end  
    end
end
 [zse,t]=reducev2(zse/n,time,400);
  [zsp,t]=reducev2(zsp/n,time,400);
   [zie,t]=reducev2(zie/n,time,400);
      [zit,t]=reducev2(zit/n,time,400);
  figure
  plot(t,zse,'b')
  hold on
  plot(t,zsp,'g')
    plot(t,zie,'r')
      plot(t,zit,'m')
    plot(t,1-zse-zsp-zie-zit,'y')
  %plot(time,1-z-zz)
% plot(time,z)
end

