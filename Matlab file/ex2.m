%==========================================================================
%   TP :            Case study: Exercse 2
%   Contact:        ezequiel.gonzalezdebada@epfl.ch
%                   matin.macktoobian@epfl.ch
%==========================================================================
classdef ex2
    %Class gathering the solutions of exercise 2. 
    methods(Static)
        %
        function [Q1,Q2] = getLQRCostFunctArrays()
            % [Q1,Q2] = getLQRCostFunctArrays(NSTATES,NINPUTS) returns the
            % matrices Q1 and Q2 defining the LQR's cost function.  
            
            Q1 = 0.5.*eye(5);
            Q1(1,1) = 0.00001;
            Q1(2,2) = 50;
%             Q1(3,3) = 5000;
            Q2 = [1, 0; 0, 0.00002];
%             Q2(2,2) = 1;
            %%- set Q1 and Q2 values. 
            %Q1 = 
            %Q2 = 
        end
        %        
        function [LQR_gain, SInf, singular_values] = getLQRGain(Phi,Gamma,Q1,Q2)
            % [LQR_K,SInf,SINGULAR_VALUES] =
            % getLQRGain(PHI,GAMMA,Q1,Q2) returns the
            % LQR control gain, the matrix SINF as well as the history of
            % the singular values of the S matrix over the iterations. The
            % method takes as inputs the PHI and GAMMA metrices describing 
            % the dicrete lineal model of the systme and the matrices Q1
            % and Q2 defining the cost function of the control problem. 
            % The last input is a string the specifies the method to be used. 
            %
            %TODO: implement algorithm in Section 7.1.7. 
            %
            % Hint: regarding the number of iterations for the calculation
            % of the control gain following Riccati approach, you should
            % find (by trial and error) a number of iterations that is
            % sufficiently large so that the S matrix converges.
            % This condition can be observed by keeping track
            % of the elements of the matrix or any other characteristics.
            % In this case, we are to use the singular values of S as it is
            % updated over the iterations.
            
%             LQR_gain = [];
%             SInf = [];
%             singular_values = [];
            
            %%- Implement the algorithm to calculate the LQR control
            %%gain using riccati equation. 

            %%- Initialize S
            S = Q1; 

            %%- number of interations
            %n_iterations = 300;
            iter = 1;
            tolerance = 1e-3;
            tol_iter = 0;

            % Vector of singular values (to control the convergence of the method
            singular_values = [];

            %Implementation of Ricatti method 
            while tol_iter <= 10

                %%- Ricatti Equation 
                R = Q2 + Gamma'*S*Gamma;
                M = S - S*Gamma*pinv(R)*Gamma'*S;
                S = Phi'*M*Phi + Q1;

                %%- store singular values (use svd function) 
                singular_values = [singular_values; svd(S)'];
                if iter > 1
                    if max(abs((singular_values(iter,:)-singular_values(iter-1,:))./singular_values(iter,:))) < tolerance
                        tol_iter = tol_iter + 1;
                    end
                end
                iter = iter +1;
            end

            %%-calculate the LQR control gain. 
            LQR_gain = pinv(R)*Gamma'*S*Phi; 
            SInf=S;


        end
        % OBSERVER-related methods        
        function Cprim = alternativeSystemsOutputEquation(C)
            % Cprim, Dprim = alternativeSystemsOutputEquation(C)
            % takes as input the C array you initially describing the 
            % system's output equation and returns the Cprime array 
            % describing the alternative output equation we will consider 
            % to implement our observer. 
            Cprim = [[eye(2);zeros(2,2)], zeros(4,1),[zeros(2,2);eye(2)]] *C;
        end
        %
        function n_states_not_observable = checkObservability(Phi,Cprime)
            % [N_STATES_NOT_OBSERVABLE] = checkObservability(PHI,Cprime)
            % checks whether the system described by Phi, 
            % Cprime is observable and returns the number of states that
            % cannot be observed. 
            %
            % Remember that the number of the non-onbservable states 
            % is calculated from the number of states and the rank
            % of the observability matrix. 
            %
            % HINT: To calculate the observability matrix, you can use the
            % matlab command 'obsv'. 
            

            %%- calculate the observability matrix
            Observability_matrix = obsv(Phi,Cprime);
            
            %%- Compute the number of non-observable states of the system
            n_states_not_observable = size(Phi,1)-rank(Observability_matrix);
        end
        %
        function [L, selected_poles] = getObserverGain(Phi, Gamma, lqr_K, Cprime)
            % [L, SELECTED_POLES] =
            % getObserverGain(obj,PHI,GAMMA,LQR_K,Cprime) 
            % returns the observer matrix L as well as the poles
            % SELECTED_POLES characterizing the observation closed loop
            % given the discrete linear model defined by PHI and GAMMA, the
            % LQR control gain LQR_K and the matrix Cprime. 
            %
            % 
            % Hint1: as a rule of thumb, the poles of the observer are 
            % calculated as some percentage of the poles that characterize 
            % the close loop dynamic. 
            %
            % Hint2: remember that the eigenvalues of the control closed
            % loop are given by eig(\Phi - \Gamma * lqr_K)
            %
            % Hint3: to calculate the observer gain, use the matlab command
            % 'place'. Read the help information and note that the
            % function is meant to calculate control gains. That
            % is it returns gains to impose certain poles assuming the
            % closed loop dynamics is given by (A-B*K). However, we are
            % using it to place the poles of the observer. Check and compare
            % the observer close loop, build a paralelism with the control
            % case, and introduce the needed modifications in the matrices
            % given to 'places'. 
            %
            % HINT4: Think about what happens when you transpose the closed
            % loop dynamic of the observer.
            
            
            %%- calculate the place where the close loop observation poles
            %%are wanted to be. 
            
            control_close_loop_poles = eig(Phi-Gamma*lqr_K);
%             percentage = 0.999;
            percentage = 0.1;
            selected_poles = percentage .* control_close_loop_poles;

            
            %%- calculate the observer gain. 
            L = place(Phi',Cprime',selected_poles)';
        end        
        %
        function x0Obs = getObserverInitialState(x_bar)
            % X0OBS = getObserverInitialState(X_BAR)
            % returns the initial guess X0OBS our observer will assign to 
            % the states it has to estimate. 
            %
            % The outputs should be column vectors, and does not need to be 
            % equal to the real initial state of the system. 
            %
            % Remember that by definition \tilde{x} = x - \overline{x} and
            % that the observer is based on the linear model. 
            x0Obs = zeros(size(x_bar,2)-1,1);
%             x0Obs(3) = -0.5;
            
%             x0Obs = {x0Obs_ex1'};
        end
        %
        function [AObs, BObs, CObs, DObs] = getObserverImplementationArrays(Phi, Gam, L, Cprime)
            % [AObs, BObs, CObs, DObs] = getObserverImplementationArrays(Phi, Gam, L, Cprime)
            % returns matrices AObs, BObs, CObs, DObs, needed to
            % implement the observer using state-space equations 
            %
            % x(k+1) = AObs x(k) + BObs u(k)
            % y(k) = CObs x(k) + DObs u(k)
            %
            % as it is illustrated in the solution of the exercises 6.6.1 -
            % 3.
            % Note that the xs and us used above are not the state and
            % input of our system, but rather the internal state and inputs
            % of the state-space block we will use to implement the
            % observer. 
            %
            %
            [N_y, N_x] = size(Cprime);
            N_u = size(Gam,2);
            AObs = Phi - L*Cprime;
            BObs = [L, Gam];
            CObs = eye(N_x);
            DObs = zeros(N_x,N_y+N_u);
            
            %%- set up arrays AObs, BObs, CObs, and DObs with the appropriate size. 
        end
    end
end

