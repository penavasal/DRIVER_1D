
function SOLVER(MAT_POINT)
    
    tic;

    global SOLVER VARIABLE GEOMETRY MATERIAL
     
    %----------------------------------------------------------------------
    % Initial state, matrixes and load
    %----------------------------------------------------------------------    
    [STEP,MAT_POINT,Int_var,Mat_state,GLOBAL]=init(MAT_POINT);
    
    %perturbated(MAT_POINT,Mat_state,Disp_field,Int_var,stiff_mtx,STEP);
       
    %----------------------------------------------------------------------
    % Solver for each Block
    %---------------------------------------------------------------------- 
    while STEP.BLCK<SOLVER.BLOCKS+1
        
        if STEP.BLCK>1
            if SOLVER.Output(STEP.BLCK)~=SOLVER.Output(STEP.BLCK-1)  || ...
                    isfile(SOLVER.Output(STEP.BLCK))==0
                save(SOLVER.Output(STEP.BLCK), 'GEOMETRY', 'VARIABLE', 'SOLVER', 'MATERIAL');
            end
        else
            save(SOLVER.Output(STEP.BLCK), 'GEOMETRY', 'VARIABLE', 'SOLVER','MATERIAL');
        end

        [STEP,MAT_POINT,GLOBAL]=MP_solver(STEP,MAT_POINT,Int_var,...
            Mat_state,GLOBAL);

        STEP.BLCK=STEP.BLCK+1;
        
        if STEP.BLCK>SOLVER.BLOCKS
            break;
        else
            % Update
            [Mat_state,Int_var]=VECTORS.Update_ini(...
                STEP,GLOBAL,Mat_state,Int_var,MAT_POINT);
        end
    end
    
    tfin=toc;

    disp(tfin);

end

function [STEP,MAT_POINT,GLOBAL]=MP_solver(...
            STEP,MAT_POINT,Int_var,Mat_state,GLOBAL)

    %--------------------------------------------------------------------------
    % Initialize variables
    %--------------------------------------------------------------------------
    global SOLVER
    
    BLCK=STEP.BLCK;
    STEP.dt = STEP.dt*SOLVER.time_factor(BLCK);
    STEP.t  = STEP.t + STEP.dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while STEP.t < SOLVER.Time_final(BLCK) 
        STEP.ste=STEP.ste+1;
        
        % 1. Problem
        [Mat_state,MAT_POINT]=PR.choice(Mat_state,Int_var,MAT_POINT,STEP);

        % 2. Constitutive &/O Stiffness_mat
        [Int_var,Mat_state,STEP]=...
            Constitutive.update(3,STEP,Int_var,Mat_state,MAT_POINT);
        
        % 3. Storage
        if rem(STEP.ste,SOLVER.SAVE_I)==0
            
            STEP.ste_p=STEP.ste_p+1;
            fprintf('ste_p %i \n',STEP.ste_p);
            
            GLOBAL=VECTORS.Store(GLOBAL,STEP,Mat_state,Int_var,MAT_POINT);
        end

        % 4. Save info
        if ((rem(STEP.ste/SOLVER.SAVE_I,SOLVER.SAVE_F)==0) || (STEP.FAIL==1))
            save(SOLVER.Output(BLCK),'MAT_POINT','GLOBAL','-append')
            if STEP.FAIL==1
                stop
            end  
        end
        
        % 5. Update
        [Mat_state,Int_var]=VECTORS.Update(Mat_state,Int_var);
            
        STEP.dt = STEP.dt*SOLVER.time_factor(BLCK);
        STEP.t  = STEP.t + STEP.dt;
    end
    
    % SAVING LAST VALUES
    STEP.ste_p=STEP.ste_p+1;
    fprintf('ste_p %i \n',STEP.ste_p);

    GLOBAL=VECTORS.Store(GLOBAL,STEP,Mat_state,Int_var,MAT_POINT);
    save(SOLVER.Output(BLCK),'MAT_POINT','GLOBAL','-append')

end

