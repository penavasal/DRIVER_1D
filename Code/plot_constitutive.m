
function plot_constitutive(e,str)

time_d=0.001;

load(str,'-mat','GLOBAL','GEOMETRY','MATERIAL','SOLVER');

BLCK=1;

mati=GEOMETRY.material(e);
MODEL=MATERIAL(BLCK).MODEL(mati,1);
MAT=MATERIAL(BLCK).MAT;

if MODEL>=2
    if ~isempty(MAT{19,mati})
        M = MAT{19,mati};
    else
        M=0;
    end
end

if ~isempty(MAT{16,mati})
    n0=MAT{16,mati};
    %n=1-(1-n0)/MAT_POINT{1}(e).J;
else
    n0=0;
end

e_0=n0/(1-n0);

Ps=GLOBAL.Ps(e,:);
P0=Ps(1);

if P0<500 % To kPa
    mult=1000;
else
    mult=0.001;
end
Ps=Ps*mult;
Sy_tot=GLOBAL.Sy(e,:)*mult;
Qs=GLOBAL.Qs(e,:)*mult;
Pw=GLOBAL.pw(e,:)*mult;

Es=GLOBAL.Es;
Es_p=GLOBAL.Es_p;
void_index=GLOBAL.J(e,:)*(1+e_0)-1;


% Time
ste_p=GLOBAL.ste_p-1;

total=100;
if total>ste_p
    each=1;
else
    each=round(ste_p/total);
end

%Ellipse

if MODEL>=3 && MODEL<4
    Pc_max=max(-Sy_tot(e,1:ste_p));
else
    Pc_max=max(-Ps(e,1:ste_p))*1.2;
end
b_max=M*Pc_max/2;
lim=1.2*Pc_max;

figure;

hold on
h2=plot(0,0);
h3=plot(0,0);
for i=1:each:ste_p
    
    if  MODEL>=2 && MODEL<3
        delete(h3)
        subplot(2,2,2)
        hold on
        C=-Sy_tot(i);
        h3=plot(linspace(0,lim,5),M*linspace(0,lim,5)+C,'k');
    end
    
    subplot(2,2,1)
    plot((-Es(e*4-2,1:i)-Es_p(e*4-2,1:i))*100,Qs(1:i))
    axis([0 inf 0 max(max(Qs(1:ste_p))*1.1,max(b_max,lim*M+C))])
    xlabel('\epsilon')
    ylabel('Q [kPa]')
    
    subplot(2,2,2)
    axis([0 Pc_max 0 max(max(Qs(1:ste_p))*1.1,max(b_max,lim*M+C))])
    xlabel('P [kPa]')
    ylabel('Q [kPa]')

    if MODEL>=3 && MODEL<4
    delete(h2)
    Pc=-Sy_tot(i);
    b=M*Pc/2;
    x0=Pc/2;
    t=0:0.01:pi;
    x=x0+Pc/2*cos(t);
    y=b*sin(t);
    
    h2=plot(x,y,'b');
    end
    
    plot(-Ps(1:i),Qs(1:i),'r')
    
    if SOLVER.UW
    subplot(2,2,3)
    plot((-Es(e*4-2,1:i)-Es_p(e*4-2,1:i))*100,Pw(1:i))
    axis([0 inf -inf max(Pw(1:ste_p))*1.1])
    xlabel('\epsilon')
    ylabel('P_w [kPa]')
    end

    subplot(2,2,4)
    plot(-Ps(2:i),1+void_index(2:i))
    axis([0 Pc_max min(void_index(2:ste_p))+0.5 max(void_index(2:ste_p))+1.5])
    xlabel('P [kPa]')
    ylabel('1+e')
    
    drawnow;
    pause(time_d)
end
hold off

end