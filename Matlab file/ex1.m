%==========================================================================
%   TP :            Case study: Exercse 1
%   Contact:        ezequiel.gonzalezdebada@epfl.ch
%==========================================================================
classdef ex1
    %Class gathering the solutions of exercise 1. 
    methods (Static)
        %
        function varargout = getSystemParameters
            % PARAMETERS = getSystemParameters() returns a 5-elements
            % column vector containing the value of the system parameters 
            % and the linearization point. Specifically it should contain 
            % (in the presented order):
            %   - k : value of curvature of the reference path.
            %   - car_length [m]: car's length.  
            %   - sigma_v : coefficient characterizing the dynamic of the
            %   actuator tracking the speed reference. 
            %   - sigma_S : coefficient characterizing the dynamic of the
            %   actuator tracking the steering wheel's reference position. 
            %   - spReg : value of the speed the vehicle should drive at. 
            k = 1e-10;
            L = 4;
            sigma_v = 1;
            sigma_S = 5;
            spReg = 5;
            varargout = {[k; L; sigma_v; sigma_S; spReg]};               
        end
        %
        function varargout = getLinealModelArrays(parameters)
            % [A,B,C,D] = getLinealModelArrays(PARAMETERS) returns the
            % matrices A,B,C,D characterizing the continuous-time linear
            % model of the system. The input *parameters* corresponds to 
            % the output of the method *getSystemParameters*.
            
            %%- Calculate arrays A,B,C,D
            Vref = parameters(5);
            k = parameters(1);
            A = zeros(5,5);
            A(1,2) = k*Vref; A(1,4) = 1; A(2,3) = Vref; A(3,2) = -k^2*Vref; 
            x5_bar = (16*atan(parameters(1)*parameters(2)));
            A(3,5) = (sec(x5_bar/16))^2*Vref/parameters(2)/16;
            A(4,4) = - parameters(3); A(5,5) = - parameters(4);
            B = zeros(5,2);
            B(4,1) = parameters(3); B(5,2) = parameters(4);
            C = eye(5);
            D = zeros(5,2);
            
            varargout = {A,B,C,D};
        end        
        %
        function varargout = getDiscreteLinearModel(A,B,C,D,sampling_time,method)
            % [PHI,GAM] =
            % getDiscreteLinearModel(A, B, C, D, SAMPLING_TIME, METHOD)
            % returns the PHI and GAMMA matrices characterizing the
            % discrete-time linear model of the system given the matrices
            % A,B,C,D of the continuous-time linear model and the desired 
            % SAMPLING_TIME.
            %
            % Additionally, the input METHOD is a string
            % indicating the method that should be used to calculate 
            % the matrices PHI and GAMMA. It can take values 
            % - Euler : Euler approximation as discretization method. 
            % - c2d : use the matlab command c2d. 
            
            switch method
                case 'Euler'
                    % Calculate the discrete-time linear model using the 
                    % Euler approximation.

                    Phi = sampling_time.*A+eye(5);
                    Gamma = sampling_time.*B;    
                    
                case 'c2d'
                    %%- Build continuous representation of the system with 'ss'
                    Mc = ss(A,B,C,D);
                    
                    %%- Calculate the discrete-time linear model of the system 
                    % using the command 'c2d'
                    
                    Md = c2d(Mc,sampling_time); 

                    %%- Extract from Md, the Phi and Gamma matrices. 

                    Phi = Md.A;
                    Gamma = Md.B;

                    %%-set up output of the function 
                case 'Psi'
                    N = 100;
                    Psi = eye(5)+A.*(sampling_time/(N+1));
                    for n = N:-1:2
                        Psi = eye(5) + A*Psi.*(sampling_time/n);
                    end
                    Phi = eye(5) + A.*sampling_time*Psi;
                    Gamma = Psi.*sampling_time*B;
            end
            
            varargout = {Phi, Gamma};
        end                
        %
        function varargout = getWorkingTrajectory(sampling_time, simulation_time, parameters)
            % [NOMINAL_TRAJECTORY_X, NOMINAL_TRAJECTORY_U] =
            % getWorkingTrajectory(SAMPLING_TIME, SIMULTAION_TIME,
            % PARAMETERS)  
            % outputs the NOMINAL_TRAJECTORY_X and NOMINAL_TRAJECTORY_U
            % given the SAMPLING_TIME between data points, the
            % SIMULATION_TIME up to which the trajectory has to be created,
            % and the vector PARAMETERS with the value sof tha system's
            % parameters.
            %
            % The outputs NOMINAL_TRAJECTORY_X, and NOMINAL_TRAJECTORY_U must
            % be arrays [t | \bar{x}] and [t | \bar{u}] 
            % whose first column corresponds to the timespan of
            % the data point, and following columns store the information
            % of the states and inputs at the corresponding time.
            %
            % The defined output trajectories are meant to be imported in
            % Simulink with the "From Workspace" block. If any
            % additional doubt regarding how the data should be formated,
            % read the information provided in the mentioned simulink block.
            %
            % todos
            % - create time vector. 
            % - create the nominal states trajectory output
            % - create the control inputs nominal trajectory output
            
            %%- create time vector
            time_vector = (0:sampling_time:simulation_time)';
            
            %%-create nominal state trajectory. 
            x1_bar = parameters(5).*time_vector;
            x2_bar = 0.*time_vector;
            x3_bar = 0.*time_vector;
%             if x4(0) = V_ref
            x4_bar = parameters(5).*ones(length(time_vector),1);
%             if x4(0) = 0
%             x4_bar = (parameters(5)*parameters(3)).*time_vector./...
%                 (parameters(3).*time_vector+ones(length(time_vector),1));
            x5_bar = (16*atan(parameters(1)*parameters(2))).*ones(length(time_vector),1);
            bar_x = [x1_bar, x2_bar, x3_bar, x4_bar, x5_bar];
            nominal_trajectory_x = [time_vector, bar_x ];
            
            %%-create nominal control input trajectory. 
            u1_bar = parameters(5).*ones(length(time_vector),1);
            u2_bar = x5_bar;
            bar_u = [u1_bar, u2_bar];
            nominal_trajectory_u = [time_vector, bar_u];
            
            varargout = {nominal_trajectory_x, nominal_trajectory_u};
        end
        %
        function varargout = getInitialState(nominal_trajectory_x)
            %[X0, X0TILDE] = getInitialState(NOMINAL_TRAJECTORY_X)
            % returns the initial state X0 of the system and the
            % initial state X0TILDE of the linear models given the 
            % information on the exercise handout and the
            % NOMINAL_TRAJECTORY_X.
            %
            % The outputs should be column vectors. 
            %
            % Remember that by definition \tilde{x} = x - \overline{x}.
            
            
            %%- define the value of x0 for experiment 1
            x0_experiment_1 = [0; 0; 0; nominal_trajectory_x(1,5); nominal_trajectory_x(1,6)];
            x0_exp_3 = [0; 0; 0; 0; 0];
            x0_exp_4 = [0; 0; pi/4; nominal_trajectory_x(1,5); nominal_trajectory_x(1,6)];
            x0_exp_5 = [0; 0; pi/4; nominal_trajectory_x(1,5); nominal_trajectory_x(1,6)];
            
            %%- define the value of x0Tilde for experiment 1
            x0Tilde_experiment_1 = x0_experiment_1 - nominal_trajectory_x(1,2:6)';
            x0Tilde_exp_3 = x0_exp_3 - nominal_trajectory_x(1,2:6)';
            x0Tilde_exp_4 = x0_exp_4 - nominal_trajectory_x(1,2:6)';
            x0Tilde_exp_5 = x0_exp_5 - nominal_trajectory_x(1,2:6)';
            
            %including the different values for different experiments as a
            %cell
%             x0 = {x0_experiment_1,x0_experiment_1,x0_exp_3,x0_exp_4,x0_exp_5};
%             x0Tilde = {x0Tilde_experiment_1,x0Tilde_experiment_1, x0Tilde_exp_3,x0Tilde_exp_4,x0Tilde_exp_5};
            
            x0 = {x0_experiment_1};
            x0Tilde = {x0Tilde_experiment_1};
            
            %set outputs of the function 
            varargout = {x0,x0Tilde};
        end
        %
        function varargout = getOpenLoopInputSignal(sampling_time, simulation_time)
            %[INPUT_CONTROL_ACTIONS_OPEN_LOOP] = getOpenLoopInputSignal(SAMPLING_TIME, SIMULATION_TIME)
            % outputs an input sequence to be applied in open loop to the 
            % models. The desired SAMPLING_TIME between data points as
            % well as the SIMULTION_TIME are provided. 
            %
            % As in the case of GETWORKINGTRAJECTORY function, the outputs
            % are meant to be used in Simulink with "From Workspace"
            % importing modules. If any additional doubt regarding how the
            % data should be structured, read the information provuded by
            % the mentioned simulink block. 
            %
            % todo:
            % - Declare an appropriate timespan vector. 
            % - Create the input_control_actions_open_loop array with the
            % sequence of control inputs to be applied in open loop. 
            %
            %
            % Notice: alternatively, this function can output a cell with
            % several arrays showing different control sequences to be
            % applied. This would make the provided script to run the
            % simulink mdel as many times as different control sequences
            % are gathered within the cell. Meaning that several
            % experiments can be set at once. 
            %
            
            %%- Create a time vector.
            %time_vector = ;
            
            %%- set the control sequence to be applied in open loop for the
            %%1st experiment. 
            time_vector = (0:sampling_time:simulation_time)';
            r1 = 0.3; r2 = pi/30;
            r3 = 5; r4 = pi/3;
            u1_open_loop = [2*r1.*rand(length(time_vector),1)+(5-r1).*ones(length(time_vector),1),...
                (2*r2).*rand(length(time_vector),1)+(16*atan(4e-10)-r2).*ones(length(time_vector),1)];
            u2_open_loop = [2*r3.*rand(length(time_vector),1)+(5-r3).*ones(length(time_vector),1),...
                (2*r4).*rand(length(time_vector),1)+(16*atan(4e-10)-r4).*ones(length(time_vector),1)];
            u4_open_loop = [3.*ones(length(time_vector),1),...
                (2*r2).*rand(length(time_vector),1)+(16*atan(4e-10)-r2).*ones(length(time_vector),1)];
            u5_open_loop = [7.*ones(length(time_vector),1),...
                (2*r2).*rand(length(time_vector),1)+(16*atan(4e-10)-r2).*ones(length(time_vector),1)];
            uOpenLoop_experiment_1 = [time_vector, u1_open_loop];
            uOpenLoop_experiment_2 = [time_vector, u2_open_loop];
            uOpenLoop_experiment_4 = [time_vector, u4_open_loop];
            uOpenLoop_experiment_5 = [time_vector, u5_open_loop];
            
            %Include different values for the different experiments in a
            %cell.
%             input_control_actions_open_loop = {uOpenLoop_experiment_1,uOpenLoop_experiment_2,...
%                 uOpenLoop_experiment_2,uOpenLoop_experiment_4, uOpenLoop_experiment_5};
            input_control_actions_open_loop = {uOpenLoop_experiment_1};

            
            %set output of the function
            varargout = {input_control_actions_open_loop};
        end
        %
    end
    
end

