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
                    [Ps,Qs]=VECTORS.invar(Mat_state.Sigma(:,1),1); 
                    [xi1,xi2]=calcula_xi2(Ps,Qs,Int_var.dgamma(1,1));  
                    S=def/steps;

                    V=-xi1 + xi2*S;

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
        
        function [xi1,xi2]=calcula_xi2(p,q,dgamma,a1)

            global MATERIAL GEOMETRY  alphag khar ghar Mg

            MODEL=MATERIAL(1).MODEL;
            Mat=GEOMETRY.material;
            
            e=1;
            
            if MODEL(Mat(e),1)<1
            
            elseif MODEL(Mat(e),1)>=4 && MODEL(Mat(e),1)<5
                
                [xi1,xi2]=PR.xi_PZ(p,q,dgamma,a1);
                
            end
        end

        
        
        function [xi1,xi2]=xi_PZ(p,q,dgamma,a1)
            
            global MATERIAL GEOMETRY
            
            Mat=GEOMETRY.material;
            MAT=MATERIAL(1).MAT;
             
            e=1;

            eta=abs(q/p);

            khar  = MAT(29,Mat(e));%khar;
            ghar  = MAT(4,Mat(e));%ghar;
            Mg    = MAT(32,Mat(e));%Mg;
            alphag= MAT(34,Mat(e));%alphag;
            
            K=khar*p;
            J=khar*q;
            G=ghar*p+(khar*(q^2)/(3*p));

            signq=sign(q);

            % Vectors
            [ng,~]=PR.build_vector(alphag,Mg,eta,q,signq);

            if (a1==-1 && q>0 ) || (a1==1 && q<0 )
                ng(1)=-abs(ng(1));
            else
                ng(1);
            end

            xi2=3*(G-J)/(3*K-J);

            xi1=dgamma * (ng(1) - xi2*ng(2));
        
        end
        
        function [n,d]=build_vector(alpha,Mf,eta,q,signq)

    %         if eta>0.1 
                if signq>=0
                    ns  = 1;
                else
                    ns  = -1;
                end
    %         else
    %             ns=0;
    %         end
            nv  = (1+alpha)*(Mf-eta); 
            d   = nv;

            det = sqrt(ns*ns+nv*nv);

            nt  = -1/2*q*Mf;

            n   = [nv ; ns; nt]/det;

        end
       
        
    end
 end
        