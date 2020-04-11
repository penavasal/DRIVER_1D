
function drained_PZ(e,str)

    load(str,'-mat','ste_p','GLOBAL');

    Ps=GLOBAL.Ps;
    Qs=GLOBAL.Qs;
    Es=GLOBAL.Es;
    
    i=ste_p;
    
    edev=zeros(i,1);
    evol=zeros(i,1);
    for j=1:i
        [edev(j),evol(j)]=invar(Es(:,j),e,4);
    end

    figure;
    subplot(2,2,1)    
    plot(edev(1:i)*100,Qs(e,1:i)./Ps(e,1:i))
    xlabel('\epsilon_s %')
    ylabel('ETA')

    subplot(2,2,2)
    plot(Ps(e,1:i)*1000,Qs(e,1:i)*1000,'r')
    xlabel('P [kPa]')
    ylabel('Q [kPa]')
    
    subplot(2,2,3)
    plot(edev(1:i)*100,evol(1:i)*100)
    xlabel('\epsilon_s %')
    ylabel('\epsilon_v %')

    subplot(2,2,4)
    plot(Qs(e,1:i)./Ps(e,1:i),evol(1:i)*100)
    xlabel('ETA')
    ylabel('\epsilon_v %')
    
end

function [Q2,P2]=invar(Ss,e,dims)

    ss=zeros(dims,1);
    for i=1:dims
        ss(dims+1-i)=Ss(e*dims+1-i);
    end

    Sc=e2E(ss,2);

    P2=(Sc(1,1)+Sc(2,2)+Sc(3,3));
    s=Sc-P2/3*eye(3);
    Q2= s(1,1)^2 + s(2,2)^2 + s(3,3)^2 +...
       s(2,1)^2 + s(1,2)^2 + ...
       s(3,1)^2 + s(1,3)^2 + ...
       s(2,3)^2 + s(3,2)^2;
    rj2=Q2*0.5;
    Q2= sqrt(2/3*Q2);
    
                %Sign
            rj3=0;
            rj3 = rj3 + 3*s(2,1)*s(2,1)*(s(1,1)+s(2,2));
            for i=1:3
                rj3 = rj3 + s(i,i)*s(i,i)*s(i,i);
            end
            rj3=rj3/3;

            rj23=sqrt(rj2)^3;
            if rj23<1.0e-15
                sint3=0;
            else
                sint3 = -3 * sqrt(3) * rj3/2/rj23;
            end
            if sint3<-1
                sint3=-1;
            elseif sint3>1
                sint3=1;
            end
            theta = 1/3*asin(sint3);
            
            if sint3<0
                Q2=-Q2;
            end
            
            P2=-P2;
end

% STRESS - STRAIN vector to tensor
function [E]=e2E(e,sp)

    if sp==1
        E=e(1);
    else
        E=zeros(3,3);
        %Build matrix
        if sp==2
            E(1,1)=e(1);
            E(2,2)=e(2);
            E(1,2)=e(4);
            E(2,1)=e(4);
            E(3,3)=e(3);
        else
            E(1,1)=e(1);
            E(2,2)=e(2);
            E(1,2)=e(4);
            E(2,1)=e(4);
            E(1,3)=e(5);
            E(3,1)=e(5);
            E(2,3)=e(6);
            E(3,2)=e(6);
            E(3,3)=e(3);
        end
    end
end
