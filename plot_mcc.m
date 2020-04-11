
function plot_mcc(e,str,M,e_0,time_d,pplot,each)

load(str,'-mat','GLOBAL');

Sy_tot=GLOBAL.Sy;
Ps=GLOBAL.Ps;
Qs=GLOBAL.Qs;
Es=GLOBAL.Es;
Es_p=GLOBAL.Es_p;
%Pw=GLOBAL.pw;
Pw=Ps*0;
void_index=GLOBAL.J*(1+e_0)-1;

ste_p=GLOBAL.ste_p

if strcmp(pplot,'FULL')
    ini=1;
elseif strcmp(pplot,'END')
    ini=ste_p;
end

%Ellipse

Pc_max=max(-Sy_tot(e,:))*1000;
lim=0.6*Pc_max;
b_max=M*Pc_max/2;

figure;

subplot(2,2,2)
hold on
h2=plot(0,0);
plot(0:lim,M*(0:lim),'k');
for i=ini:each:ste_p
    
    subplot(2,2,1)
    plot(-Es(e*4-2,1:i)-Es_p(e*4-2,1:i),Qs(e,1:i)*1000)
    axis([0 inf 0 max(max(Qs(e,:))*1.1*1000,b_max)])
    xlabel('\epsilon')
    ylabel('Q [kPa]')
    
    subplot(2,2,2)
    axis([0 Pc_max 0 max(max(Qs(e,:))*1.1*1000,b_max)])
    xlabel('P [kPa]')
    ylabel('Q [kPa]')

    delete(h2)
    Pc=-Sy_tot(e,i)*1000;
    b=M*Pc/2;
    x0=Pc/2;
    t=0:0.01:pi;
    x=x0+Pc/2*cos(t);
    y=b*sin(t);
    
    h2=plot(x,y,'b');
    plot(-Ps(e,1:i)*1000,Qs(e,1:i)*1000,'r')
    
    subplot(2,2,3)
    plot(-Es(2,1:i)-Es_p(2,1:i),Pw(1,1:i)*1000)
    axis([0 inf -inf max(Pw(e,:))*1.1*1000])
    xlabel('\epsilon')
    ylabel('P_w [kPa]')

    subplot(2,2,4)
    plot(-Ps(e,2:i)*1000,1+void_index(e,2:i))
    axis([0 Pc_max min(void_index(e,2:ste_p))+0.5 max(void_index(e,2:ste_p))+1.5])
    xlabel('P [kPa]')
    ylabel('1+e')
    
    drawnow;
    pause(time_d)
end
hold off

end