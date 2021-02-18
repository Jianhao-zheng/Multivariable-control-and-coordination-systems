%==========================================================================
%   TP :            Case study: Exercse 1
%   Contact:        ezequiel.gonzalezdebada@epfl.ch
%                   matin.macktoobian@epfl.ch
%==========================================================================

%% Clear the workspace and close figures
%==========================================================================
clear all, close all, clc;

%% INITIALIZATION 
%==========================================================================
    animate = false;
    sampling_time = 0.1;
    simulation_time = 10;
    ut = utilities;
    sol = ut.load_solution_class(1);

%% Get system parameters
%==========================================================================
    parameters = sol.getSystemParameters;
    if ut.stopCondition(parameters, 'getSystemParameters'), return; end
    
%% Construct the path to be tracked
%==========================================================================    
    curvature = inf;
    if ~isempty(parameters)
        curvature = parameters(1);
    end
    path_to_track = ut.genCircularPath(curvature);

    %% Get arrays defining the linear model (calling your solutions). 
    %==========================================================================
    [A,B,C,D] = sol.getLinealModelArrays(parameters);
    if ut.stopCondition(A, 'getLinealModelArrays'), return; end

    % Calculate discrete models
    [Phi_euler1, Gam_euler] = sol.getDiscreteLinearModel(A,B,C,D,sampling_time,'Euler');
    [Phi,Gam] = sol.getDiscreteLinearModel(A,B,C,D,sampling_time,'c2d');
    [Phi_psi,Gam_psi] = sol.getDiscreteLinearModel(A,B,C,D,sampling_time,'Psi');
%     Phi = Phi_euler1;
%     Gam = Gam_euler;
    % Printing discrete model matrices in command window
    if ut.stopCondition(Phi_euler1, 'getDiscreteLinearModel'), return; end

    %Print Matrices
    ut.printLabeledAB(A, B, 'Continuous system','description','[A | B]');
    ut.printLabeledAB(Phi_euler1, Gam_euler, 'Discrete description obtained applying Euler approximation');
%     ut.printLabeledAB(Phi_psi, Gam_psi, 'Discrete description obtained programmatically');
    ut.printLabeledAB(Phi, Gam, 'Discrete description obtained using Matlab functions');

    % Comparing Euler to Matlab
    ut.printRelativeError(Phi, Phi_euler1, Gam, Gam_euler, 'Comparing Euler vs Matlab')
%     ut.printRelativeError(Phi, Phi_psi, Gam, Gam_psi, 'Comparing Programmatically vs Matlab')
    
    % time vector (required to generate nominal trajectories)
    nominal_trajectory_sampling_time = sampling_time/100;

    % Generate nominal trajectories for simulation
    [nominal_trajectory_x, nominal_trajectory_u] = sol.getWorkingTrajectory(nominal_trajectory_sampling_time, simulation_time, parameters);
    if ut.stopCondition(nominal_trajectory_x, 'getWorkingTrajectory'), return; end

%% Simulations
%==========================================================================
    % Considering the nominal trajectory, get initial state. 
    [x0_experiments, x0Tild_experiments] = sol.getInitialState(nominal_trajectory_x);
    if ut.stopCondition(x0_experiments, 'getInitialState'), return; end

    % Generate open loop control reference 
    [input_sequence_exps] = sol.getOpenLoopInputSignal(nominal_trajectory_sampling_time,simulation_time);
    if ut.stopCondition(input_sequence_exps, 'getOpenLoopControlSignal'), return; end

    for exp_id = 1:length(input_sequence_exps)
        fprintf('\nExperiment %02d out of %02d\n', exp_id, length(input_sequence_exps));

        x0 = x0_experiments{exp_id};
        x0_tild = x0Tild_experiments{exp_id};
        input_sequence = input_sequence_exps{exp_id};

        %run simulate 
        sim(ut.complete_simulink_model_name('open_loop_experiments'));

        % plot results
        ax_cells = ut.plotResults(time, discrete_time, continuous_non_linear_model_global_position, open_loop_u, continuous_non_linear_model_y, path_to_track, parameters, 'figure_name', sprintf('open Loop - exp%02d: state',exp_id),'plot_label','continuous-time non-linear model');
        ut.plotResults(time, discrete_time, continuous_linear_model_global_position, open_loop_u, continuous_linear_model_y, path_to_track, parameters, 'plot_label','continuous-time model','ax_cells',ax_cells,'animate',animate);
        ut.plotResults(time, discrete_time, discrete_linear_model_global_position, open_loop_u, discrete_linear_model_y, path_to_track, parameters, 'plot_label','discrete-time linear model','ax_cells',ax_cells,'animate',false);
    end