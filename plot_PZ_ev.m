
function plot_PZ_ev(e,str)

    load(str,'-mat','ste_p','GLOBAL');

    Es=GLOBAL.Es;
    Esp=GLOBAL.Es_p;
    
    i=ste_p;
    
     edev=zeros(i,1);
     evol=zeros(i,1);
     for j=1:i
         [edev(j),evol(j)]=invar2(Es(:,j)+Esp(:,j),e,4);
     end
    

    figure;
    
    plot(edev(1:i)*100,-evol(1:i)*100/3)
    xlabel('\epsilon_s %')
    ylabel('\epsilon_v %')
           
end

function [Q2,P2]=invar2(Ss,e,dims)

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
