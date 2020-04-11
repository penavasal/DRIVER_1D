
function plot_PZ(e,str,e_0)

    load(str,'-mat','GLOBAL');

    Ps=GLOBAL.Ps;
    Qs=GLOBAL.Qs;
    Es=GLOBAL.Es;
    Esp=GLOBAL.Es_p;
    Pw=GLOBAL.pw;
    ste_p=GLOBAL.ste_p;

    void_index=GLOBAL.J*(1+e_0)-1;
    
    i=ste_p;
    
     edev=zeros(i,1);
     for j=1:i
         [edev(j)]=invar2(Es(:,j)+Esp(:,j),e,4);
     end
    

    figure;
        
    subplot(2,2,1)
    plot(edev(1:i)*100,Qs(e,1:i)*1000)
    xlabel('\epsilon %')
    ylabel('Q [kPa]')

    subplot(2,2,2)
    xlabel('P [kPa]')
    ylabel('Q [kPa]')
    plot(-Ps(e,1:i)*1000,Qs(e,1:i)*1000,'r')

    subplot(2,2,3)
    plot(edev(1:i)*100,Pw(1,1:i)*1000)
    xlabel('\epsilon %')
    ylabel('P_w [kPa]')

    subplot(2,2,4)
    plot(-Ps(e,2:i)*1000,1+void_index(e,2:i))
    xlabel('P [kPa]')
    ylabel('1+e')
    
end

function [Q2]=invar2(Ss,e,dims)

    ss=zeros(dims,1);
    for i=1:dims
        ss(dims+1-i)=Ss(e*dims+1-i);
    end

    Sc=e2E(ss,2);

    P2=(Sc(1,1)+Sc(2,2)+Sc(3,3))/3;
    s=Sc-P2*eye(3);
    Q2= s(1,1)^2 + s(2,2)^2 + s(3,3)^2 +...
       s(2,1)^2 + s(1,2)^2 + ...
       s(3,1)^2 + s(1,3)^2 + ...
       s(2,3)^2 + s(3,2)^2;
    Q2= sqrt(2/3*Q2);
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
