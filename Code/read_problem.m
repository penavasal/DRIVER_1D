
function MAT_POINT=read_problem(str)

% File: read_problem
%   Read some important problem parameters from problem.txt
%
% Date:
%   Version 3.0   26.11.2019

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some possible parameters if they are not read
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global VARIABLE SOLVER PROBLEM
     
    % GEOMETRY    
        
    SOLVER.AXI=0;
    
    SOLVER.UW=0;
       
    SOLVER.INIT_STEP=0;
    
    SOLVER.TYPE={'' ''};
    
    SOLVER.thickness=1;
    
    SOLVER.FAIL = 0;
    
    SOLVER.BLOCKS = 1;
    
    SOLVER.SMALL = 0;
    
    SOLVER.FRAC = 0;
    
    SOLVER.PHASES = {'U' '1'};
    
    % VARIABLES
    VARIABLE.g=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read problem.txt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen(str, 'rt'); % opción rt para abrir en modo texto
    formato = '%s %s %s %s %s %s %s %s %s'; % formato de cada línea 
    data = textscan(fid, formato, 'HeaderLines', 1);
    a = data{1};
    % Convertir a vector numérico
    b1= cellfun(@str2num, data{2}, 'UniformOutput', false);
    b2= data{2};
    b3= data{3};
    b4= data{4};
    [l,~] = size(b1);
    b=zeros(l,1);
    for i=1:l
        if b1{i}
            b(i)=b1{i};
        else
            b(i)=0;
        end
    end
    
    
    SOLVER.BLOCKS=1;
    
    SOLVER.time_factor=1*ones(SOLVER.BLOCKS,1);
    FILES=strings(SOLVER.BLOCKS,4);
    
    PROBLEM.TYPE=strings(2,1);
    PROBLEM.Q=zeros(1,2);
    PROBLEM.DEF=0;
    PROBLEM.LOOP=0;
    
    OUT=strings(SOLVER.BLOCKS,1);
    
    SOLVER.OutputType=zeros(1,2);

    BLCK=1;
    t=0;
    while (t<l)
        t=t+1;
        s1=a{t};
        
        if s1=='*'
            if strcmp(b2{t},'TYPE_OF_PROBLEM')
                disp('Reading problem')
                continue
            else
                break
            end
        end
        
        switch s1
            case '//'
                continue
            case 'INIT_FILE'
                val = str2double(b2{t});
                if isnan(val)
                    SOLVER.INIT_file=b2{t};  % Start from a previous stage
                else
                    SOLVER.INIT_file=b(t);
                end
                continue
            case 'INIT_STEP'
                SOLVER.INIT_STEP=b(t);  %
                continue 
            case 'FRAMEWORK'
                if strcmp(b2{t},'SMALL_STRAIN')
                    SOLVER.SMALL=1;
                end
                continue
            case 'SAVE_FREQUENCY'
                SOLVER.SAVE_I=b(t); % Save info each XX steps
                continue
            case 'FILE_FREQUENCY'
                SOLVER.SAVE_F=b(t); % Save the file each XX steps
                continue
            case 'TIME_FINAL'
                SOLVER.Time_final(BLCK)=b(t);
                continue
            case 'TIME_STEP'
                SOLVER.time_step(BLCK)=b(t);  
                continue
            case 'TIME_FACTOR'
                SOLVER.time_factor(BLCK)=b(t);  
                continue
            case 'OUTPUT'
                OUT(BLCK,1)=b2{t};
                continue
            case 'MATERIAL'
                FILES(BLCK,1)=b2{t};
                continue
            case 'PROBLEM_TYPE'
                PROBLEM.TYPE(1,BLCK)=b2{t};
                continue  
            case 'DRAINED'
                if strcmp(b2{t},'YES')
                    PROBLEM.TYPE(2,BLCK)='DRAINED';
                elseif strcmp(b2{t},'NO')
                    PROBLEM.TYPE(2,BLCK)='UNDRAINED';
                else
                    error('Error, unrecognized drained parameter: %s !!\n',s1)
                end
                continue 
            case 'FINAL_DEVIATORIC_STRAIN'
                PROBLEM.DEF(1,BLCK)=b(t);
                continue 
            case 'STEPS_BY_CYCLE'
                PROBLEM.LOOP(1,BLCK)=b(t);
                continue 
            case 'MAX_Q'
                PROBLEM.Q(BLCK,1)=b(t);
                continue 
            case 'MIN_Q'
                PROBLEM.Q(BLCK,2)=b(t);
                continue 
            otherwise
                fprintf('Error, unrecognized parameter: %s !!\n',s1)
                stop
        end
    end
    
%     % BLOCKS
%     if strcmp(b2{t},'NUMBER_OF_BLOCKS')
%         t=t+1;
%         SOLVER.BLOCKS=str2double(a{t});
%         fprintf('Number of blocks: %i \n',SOLVER.BLOCKS)
%     else
%         fprintf('Error, unrecognized parameter: %s !!\n',s1)
%         stop
%     end
    
    fclose(fid); 

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add the geometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MAT_POINT=geometry(SOLVER.AXI);
    
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BLOCK problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:SOLVER.BLOCKS
        
        if strcmp(OUT(i),'')
           SOLVER.Output(i)='FILE.mat';
        else
           SOLVER.Output(i)=strcat(OUT(i),'.mat');
        end

        A=exist(SOLVER.Output(i),'file');

        j=0;
        while A==2
            j=j+1;
            SOLVER.Output(i)=strcat(OUT(i),'_',int2str(j),'.mat');
            A=exist(SOLVER.Output(i),'file');
        end
        
        %--------------------------------------------------------------------------
        % PARAMETERS
        %--------------------------------------------------------------------------
        % Read material properties
        read_material(FILES(i,1),i);
        
         %**************
         % SOLVER
         %**************
 
         % Time-related conditions
         %----------------------------------------------------------------------
         time_variables(i); 
    end
end

function time_variables(i)
    global SOLVER

    % CFL CALCULATION
    %------------------------------------------------
    dt=SOLVER.time_step(i);
    tf=SOLVER.time_factor(i);

    %ste=floor(SOLVER.Time_final(i)/tf/dt)+1;
    tb=0;
    ste=0;
    while tb<SOLVER.Time_final(i)
        ste=ste+1;
        tb=tb+dt*tf;
        dt=dt*tf;
    end


    if i==1
        SOLVER.step_final(i)=ste;
        SOLVER.Time_final(i)=tb;
    else
        SOLVER.step_final(i)=ste+SOLVER.step_final(i-1);
        SOLVER.Time_final(i)=tb+SOLVER.Time_final(i-1);
    end

    if SOLVER.SAVE_I==1
        SOLVER.dim=floor(SOLVER.step_final(i)/SOLVER.SAVE_I)+SOLVER.BLOCKS;
    else
        SOLVER.dim=floor(SOLVER.step_final(i)/SOLVER.SAVE_I)+1+SOLVER.BLOCKS;
    end

    fprintf('%i plot steps when finish BLOCK %i\n',SOLVER.dim,i);
    fprintf('Save %i times when finish BLOCK %i\n',...
        round(SOLVER.dim/SOLVER.SAVE_F),i);
end

function MAT_POINT=geometry(AXI)

    global GEOMETRY

    GEOMETRY.material(1)=1;
    GEOMETRY.mat_points=1;
    GEOMETRY.sp=2;
    GEOMETRY.df=2;
    GEOMETRY.Area=1;
    
    
    if AXI==0
        if GEOMETRY.sp==2
            GEOMETRY.b_dim=3;
            GEOMETRY.s_dim=4;
            GEOMETRY.f_dim=5;
        elseif GEOMETRY.sp==1
            GEOMETRY.b_dim=1;
            GEOMETRY.s_dim=1;
            GEOMETRY.f_dim=1;
        else
            GEOMETRY.b_dim=6;
            GEOMETRY.s_dim=6;
            GEOMETRY.f_dim=9;
        end  
    elseif AXI==1
        GEOMETRY.b_dim=4;
        GEOMETRY.s_dim=4;
        GEOMETRY.f_dim=5;
    else 
        disp('DIMENSION OF THE AXI PARAMETER??')
        stop
    end
    
    MAT_POINT{1}(1:GEOMETRY.mat_points)=...
          struct(...
                 'element', ones(1),...
                 'J', ones(1),...
                 'w', zeros(1),...
                 'xi', zeros(GEOMETRY.sp,1),...
                 'near_p',zeros(1)...
             );

end

