%==========================================================================
%   TP :            Case study: Exercse 2
%   Contact:        ezequiel.gonzalezdebada@epfl.ch
%                   matin.macktoobian@epfl.ch
%==========================================================================
clear all;close all;clc;

%% script inputs
%==========================================================================
    animate_flag = true;
    simulation_time = 100;
    sampling_time = 0.01;
    sampling_time_nominal_trajectory = sampling_time/2;
    
    
    ut = utilities;
    sol1 = ut.load_solution_class(1);
    sol = ut.load_solution_class(2);
    ignore_observer = 1;
    
%% load the linearized discrete model (inherited from the previous exercise)
%==========================================================================
    parameters = sol1.getSystemParameters;
    if ut.stopCondition(parameters, 'getSystemParameters'), return; end
    [A,B,Corig,D] = sol1.getLinealModelArrays(parameters);
    if ut.stopCondition(A, 'getLinealModelArrays'), return; end
    [Phi,Gam] = sol1.getDiscreteLinearModel(A,B,Corig,D,sampling_time,'c2d');
    if ut.stopCondition(Phi, 'getDiscreteLinearModel'), return; end
    % Loading the system linearization state and initial state 
    [x_bar, u_bar] = sol1.getWorkingTrajectory(sampling_time_nominal_trajectory, simulation_time, parameters);
    if ut.stopCondition(x_bar, 'getWorkingTrajectory'), return; end
    [x0s,x0Tildes] = sol1.getInitialState(x_bar);
    if ut.stopCondition(x0s, 'getInitialState'), return; end
    
%% Optimal Controller design
%==========================================================================
    % time vector (required to generate nominal trajectories)    
	path_options = {'circle', 'path_1', 'path_2', 'path_3'};
    chosen_path = path_options{4};
    load(chosen_path);
    
    [Q1,Q2] = sol.getLQRCostFunctArrays();
    fprintf('+ Retrieving Q1 and Q2...')
    if ut.stopCondition(Q1, 'getLQRCostFunctArrays'), return; end
    fprintf(' done\n')
    
    %According to Q1 and Q2, calculate LQR gain iteratively. 
    fprintf('+ Calculating LQR gain...')
    [lqr_gain, SInf, S_vs]  = sol.getLQRGain(Phi,Gam,Q1,Q2);
    [lqr_test, SInf_test, S_test] = dlqr(Phi,Gam,Q1,Q2);
    fprintf(' done\n')
    
    fprintf('+ Representing the convergency of the S array...')
    if ut.stopCondition(lqr_gain, 'calculateLQRGain - Riccati method'), return; end
    ut.drawSingularValEvol(S_vs);
    fprintf(' done\n')
    
%% Observer design 
%==========================================================================
    % obtain the new C, excluding some state from the output, 
    fprintf('+ Retrieving arrays defining the alternative output equation...')
    Cprim = sol.alternativeSystemsOutputEquation(Corig);
    fprintf(' done\n')
    
    % Check if the system is observable with the new outputs
    fprintf('+ Checking observability...')
    unobStates = sol.checkObservability(Phi,Cprim);
    if ut.stopCondition(unobStates, 'checkObservability'), return; end
    fprintf(' %d states cannot be observed with this configuration... ',unobStates);
    fprintf(' done\n')

    fprintf('+ Calculating observer gain...')
    [L, poles_options] = sol.getObserverGain(Phi,Gam,lqr_gain,Cprim);
    if ut.stopCondition(L, 'getObserverGain'), return; end
    ut.printComplexNumbers(poles_options, 'label', sprintf('  Observer poles'))
    fprintf('  done\n')

    % Design an observer using pole placement
    fprintf('+ Retriving observer'' initial guess...')
    x0Obs = sol.getObserverInitialState(x_bar);
    fprintf(' done\n')
    
    % Matrices to 
    fprintf('+ Building arrays to implement the observer using a state-space model...')
    [AObs, BObs, CObs,DObs] = sol.getObserverImplementationArrays(Phi, Gam, L, Cprim);
    if ut.stopCondition(Cprim, 'getCArrayConsideringMeasurableStates'), return; end
    fprintf(' done\n')

%%
fprintf('+ Running simulations...\n')
n_simulations = length(x0s);
for i =1:n_simulations 
    fprintf('++ Running simulation pairs %02d/%02d...\n', i, n_simulations)
    x0 = x0s{i};
    x0Tilde = x0Tildes{i};
    %%  Run simulation with fully measurable states. 
    %========================================================================== 
    fprintf(' ++ Simulating the control loop without considering the observer...\n')
    LQR_input_signal_ID = 1;
    C = Corig;
    %Simulate and represent results
    %
    sim(ut.complete_simulink_model_name('LQR_observer'));
    
    ax_cells = ut.plotResults(time, discrete_time, x_y_theta, u, x, path_to_track, parameters, 'figure_name', 'LQR - non linear - clean and fully measurable state','animate',animate_flag);
    
    %% Run simulation with partially measurable states. 
    %==========================================================================
    fprintf(' ++ Simulating the control loop including the observer...\n')
    LQR_input_signal_ID = 2;
    ignore_observer = 0;
    sim(ut.complete_simulink_model_name('LQR_observer'));
    noise = 0;
    ut.compareObservedVsRealStates(time, discrete_time, x_hat1, x ,x, 'figure_name', 'comparison observed states')
    ut.plotResults(time, discrete_time, x_y_theta, u, x, path_to_track, parameters, 'figure_name', 'LQR + observer', 'animate', animate_flag);
end