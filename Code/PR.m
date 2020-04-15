 classdef PR
    methods(Static)
        
        % Plastic Strains in Mat. Points to nodes
        function [Mat_state,MAT_POINT]=choice(Mat_state,Int_var,MAT_POINT,STEP)

            global PROBLEM SOLVER GEOMETRY
            
            type=PROBLEM.TYPE(1,1);
            subtype=PROBLEM.TYPE(2,1);
            
            def=PROBLEM.DEF;
            
            steps=SOLVER.step_final;
            
            fdim=GEOMETRY.f_dim;
            
            
            %Vector F, Fp to Matrixes
            f_v       = zeros(fdim,1);
            for i=1:fdim
                f_v(i,1)=Mat_state.F(i,1);
            end           
            [F0]=LIB.v2m(f_v);
            
            if strcmp(type,'SIMPLE_SHEAR')
                if strcmp(subtype,'UNDRAINED')

                    t=F0(1,2);
                    
                    t=t+def/steps;
                    F=[1 t 0;
                        0 1 0;
                        0 0 1];
                    %Storage of vector F
                    [f]=LIB.m2v(F);

                    for i=1:fdim
                        Mat_state.F(fdim+1-i,1)=f(fdim+1-i);  
                    end
                    MAT_POINT{1}(1).J=det(F);
                else
                    error('Simple shear is always undrained!!!')
                end
            elseif strcmp(type,'TRIAXIAL')
                if strcmp(subtype,'UNDRAINED')
                    t=F0(2,2)-1;
                    
                    t=t-def/steps;
                    r=sqrt(1/(1+t))-1;
                    
                    F=[1+r 0 0;
                        0 1+t 0;
                        0 0 1+r]; 
                else
                    S=-def/steps;
                    St=S*STEP.ste;
                    
                    [Ps,Qs]=VECTORS.invar(Mat_state.Sigma(:,1),1); 
                    [xi1,xi2]=PR.calcula_xi2(Ps,Qs,Int_var,St);  
                    
                    V= xi1 + xi2*S;

                    e1=V/3+S;
                    e3=V/3-S/2;

                    dF=[exp(e3) 0 0;
                        0 exp(e1) 0;
                        0 0 exp(e3)];

                    F=dF*F0;   
                end
                    
            elseif strcmp(type,'TRIAXIAL_CYCLIC')
                ste1=PROBLEM.LOOP;
                j=STEP.ste;
                if ste1~=0
                    aa=(j-rem(j,ste1))/ste1;
                    if mod(aa,2)==0
                        a1=-1;
                    else
                        a1=1;
                    end

                    t=F0(2,2)-1;

                    t=t+a1*def/ste1;
                    r=sqrt(1/(1+t))-1;

                    F=[1+r 0 0;
                        0 1+t 0;
                        0 0 1+r];
                else
                    Q1=PROBLEM.Q(1,1);
                    Q2=PROBLEM.Q(1,2);
                    if j==2
                        PROBLEM.a1=1;
                    end
                    [~,Qs]=VECTORS.invar(Mat_state.Sigma(:,1),1); 

                    if Qs>Q1 && PROBLEM.a1>0 
                        PROBLEM.a1=-PROBLEM.a1;
                    elseif Qs<Q2 && PROBLEM.a1<0 
                        PROBLEM.a1=-PROBLEM.a1;
                    end

                    S=-PROBLEM.a1*def/steps;
                    V=0;

                    e1=V/3+S;
                    e3=V/3-S/2;

                    dF=[exp(e3) 0 0;
                        0 exp(e1) 0;
                        0 0 exp(e3)];

                    F=dF*F0;
                end
            end

            %Storage of vector F
            [f]=LIB.m2v(F);

            for i=1:fdim
                Mat_state.F(fdim+1-i,1)=f(fdim+1-i);  
            end
            MAT_POINT{1}(1).J=det(F);

        end
        
        function [xi1,xi2]=calcula_xi2(p,q,Int_var,St)

            global MATERIAL GEOMETRY

            MODEL=MATERIAL(1).MODEL;
            MAT=MATERIAL(1).MAT;
            Mat=GEOMETRY.material;
            
            e=1;
            
            if MODEL(Mat(e),1)<3
                if MODEL(Mat(e),1)<3
                    D22 = 9*MAT{4,Mat(e)};
                else
                    D22 = 3*MAT{4,Mat(e)};
                end
                D11 = MAT{29,Mat(e)};
                D12 = 0;
                
            elseif MODEL(Mat(e),1)>=3 && MODEL(Mat(e),1)<4
                [D11,D22,D12]=PR.xi_MCC(p,q,Int_var.gamma,St);
            elseif MODEL(Mat(e),1)>=4 && MODEL(Mat(e),1)<5
                [D11,D22,D12]=PR.xi_PZ(p,q,Int_var.P0);
            end
            xi2=(D22+3*D12)/(3*D11+D12);
            xi1=(Int_var.epsv(e,2)-Int_var.epsv(e,3))+...
                xi2*(Int_var.gamma(e,2)-Int_var.gamma(e,3));
            
            e;
        end

        
        
        function [K,D22,J]=xi_PZ(p,q,P0)
            
            global MATERIAL GEOMETRY
            
            Mat=GEOMETRY.material;
            MAT=MATERIAL(1).MAT;
             
            e=1;
            
            khar  = -MAT{29,Mat(e)}/P0(1);%khar;
            ghar  = -MAT{4,Mat(e)}/P0(1);%ghar;
            
            K=-khar*p;
            J=khar*q;
            G=-(ghar*p+(khar*(q^2)/(3*p)));
            
            D22=3*G;
        end
        
        function [K,G,J]=xi_MCC(p,q,gamma,St)
            
            global MATERIAL GEOMETRY
            
            Mat=GEOMETRY.material;
            MAT=MATERIAL(1).MAT;
             
            e=1;
            
            mu0 = MAT{4,Mat(e)};   %mu0
            k   = MAT{22,Mat(e)};  %kappa
            
            epses=St+gamma(e,1);
            
            if gamma(e,1)~=0
                e;
            end
            
            % alfa = MAT{20,Mat(e)};  %alfa
            % OMEGA = -(epsev-epsev0)/kappa;
            
            K=-p/k;
            G=-q/epses;
            J=k*(G-3*mu0*epses);

        end
    end
 end
        