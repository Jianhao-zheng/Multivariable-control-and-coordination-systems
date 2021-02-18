%==========================================================================
%   Contact:        ezequiel.gonzalezdebada@epfl.ch
%==========================================================================
classdef utilities
    %UTILITIES Class gathering some auxiliary functions used during the
    %case study of the course "multivariable control and coordination
    %systems" 
    methods
        function ax_cells = plotResults(obj, varargin)
            % plotResults(TIME, DISCRETE_TIME, GLOBAL_TRAJECTORY,
            % CONTROL_INPUT, SYSTEM_STATE, PATH_TO_TRACK, PARAMETERS)
            % plots the simulation results given the simulation outputs
            % TIME, DISCRETE_TIME, GLOBAL_TRAJECTORY,
            % CONTROL_INPUT, SYSTEM_STATE, the PATH_TO_TRACK defined by the
            % user, and the PARAMETERS characterizing the system. 
            %
            % The functon outputs a cell AX_CELLS containing the axes handle of the
            % plots.
            %
            % Additionally, those inputs can be followed by parameter/value
            % pairs to specify additional properties of the script. 
            %
            %   - ax_cells : cell containing the axes handle of the plots.
            %   - plot_label : label of the plot which will appear on the
            %     legend.
            %   - figure_name : Name of the resulting figure. 
            %   - noisy_states : if noise is considered in the simulation, 
            %     the noisy states of the systme can optionally be ploted. 
            %   - animate : boolean flag showing whether an animation of
            %     the results is wanted to be shown or not. 
            p = inputParser;
            p.addRequired('time')
            p.addRequired('discrete_time')
            p.addRequired('global_trajectory')
            p.addRequired('control_inputs')
            p.addRequired('system_state')
            p.addRequired('path_to_track')
            p.addRequired('parameters')
            
            p.addOptional('ax_cells',[])
            p.addOptional('plot_label',[])
            p.addOptional('figure_name',[])
            p.addOptional('noisy_states',[])
            p.addOptional('animate',false)
            p.parse(varargin{:})
            args = p.Results;
            ax_cells = args.ax_cells;
            
            if isempty(ax_cells)
                ax_cells = cell(3,1);
                
                f1 = figure;
                n_us = size(args.control_inputs,2);
                ax_cells{1} = gobjects(n_us,1);
                lab = 'Control input';
                for i = 1:n_us
                    ax_cells{1}(i) = subplot(n_us,1,i,'Parent',f1,'NextPlot','add');
                    title(sprintf('%s - %d',lab,i));
                end                
                
                f2 = figure;
                n_xs = size(args.system_state,2);
                ax_cells{2} = gobjects(n_xs,1);
                lab = 'State';
                for i = 1:n_xs
                    ax_cells{2}(i) = subplot(n_xs,1,i,'Parent',f2,'NextPlot','add');
                    title(sprintf('%s - %d',lab,i));
                end                                

                f3 = figure;
                ax_cells{3} = gobjects(1,1);
                ax_cells{3} = axes('Parent',f3,'NextPlot','add');
                if ~isempty(args.figure_name)
                    f1.Name = sprintf('%s - inputs', args.figure_name);
                    f2.Name = sprintf('%s - states', args.figure_name);
                    f3.Name = sprintf('%s - trajectory', args.figure_name);
                end                
            end
            
            obj.drawColumnsSeparately(args.time, args.discrete_time, args.control_inputs,ax_cells{1},args.plot_label,'x_label','time','y_labels',{'v_{ref} [m/s]','\phi_{red} [rad]'},'stairs',true);
            obj.drawColumnsSeparately(args.time, args.discrete_time, args.system_state,ax_cells{2},args.plot_label,'x_label','time','y_labels',{'s [m]','d [m]','\theta_e [rad]','v [m/s]','\phi [rad]'},'noisy_y_data',args.noisy_states);
            axis(ax_cells{3},'equal');
            
            legend_texts = legend(ax_cells{3}(1));
            if isempty(legend_texts) || (isprop(legend_texts,'String') && isempty(legend_texts.String))
                plot(args.path_to_track.coordinates(:,1),args.path_to_track.coordinates(:,2),'Parent',ax_cells{3},'LineStyle','--');
                legend(ax_cells{3}(1), 'Reference path');
            end
            
            plot(args.global_trajectory(:,1),args.global_trajectory(:,2),'Parent',ax_cells{3},'LineWidth',2);
            if length(legend_texts.String)>1 && isempty(args.plot_label)
                legend_texts.String{2} = 'Vehicle path';
            end
            
            x_range = [min(args.global_trajectory(:,1)),max(args.global_trajectory(:,1))];
            y_range = [min(args.global_trajectory(:,2)),max(args.global_trajectory(:,2))];
            x_lim = mean(x_range) + [-.5 .5]*max( [(x_range(2)-x_range(1)), (y_range(2)-y_range(1))]);
            y_lim = mean(y_range) + [-.5 .5]*max( [(x_range(2)-x_range(1)), (y_range(2)-y_range(1))]);
            y_lim(2) = y_lim(2) + 1e-10 * (y_lim(2) == y_lim(1));
            x_lim(2) = x_lim(2) + 1e-10 * (x_lim(2) == x_lim(1));
            
            set(ax_cells{3},'XLim',x_lim,'YLim',y_lim);
            
            if ~isempty(args.plot_label)    
                legend_texts = legend(ax_cells{3}(1));
                if ~isempty(legend_texts), legen_strings = legend_texts.String(1:(end));
                else, legen_strings = {}; end
                strings = {legen_strings{:}, args.plot_label};
                valid_idx = ones(length(strings),1);
                for k = 1:length(valid_idx)
                    if strcmp(strings{k}(1:4),'data')
                        valid_idx(k) = 0;
                    end
                end
                strings = strings(valid_idx>0);
                legend(ax_cells{3}(1), strings{:});     
            end
            obj.animateResults(args.time, args.discrete_time, args.global_trajectory, args.path_to_track.coordinates, args.parameters,'label',args.plot_label,'plot',args.animate); 
        end
        %
        function selected_time = select_time_vector(~, time, discrete_time, vector, varargin)
            % select_time_vector(TIME, DISCRETE_TIME, VECTOR)
            % returns the approriate time vector (to be selected from
            % TIME and DISCRETE_TIME) whose dimension corresponds with the
            % dimension of VECTOR.             
            selected_time = time;
            if size(vector,1)~=length(time) && size(vector,1)==length(discrete_time)
                selected_time = discrete_time;
            end                      
        end
        %
        function animateResults(obj, time, discrete_time, global_trajectory, path, parameters, varargin)
            % animateResults(TIME, GLOBAL_TRAJECTORY, PATH,
            % PARAMETERS) animates the simulation results given the TIME
            % vector, the GLOBAL_TRAJECTORY described by the vehicle, the
            % PATH to be tracked by the vehicle and the PARAMETERS vector.
            %
            % Additionally, those inputs can be followed by soe
            % parameter/value pairs to specify additional properties of the
            % script.
            %
            %   - label : string to be set as figure title.
            %   - plot: boolean parameter showing wether the animation
            %   should be represented. 
            p = inputParser;
            p.addOptional('label','');
            p.addOptional('plot',true);
            p.parse(varargin{:});
            args = p.Results;
            
            if args.plot
                selected_time = obj.select_time_vector(time, discrete_time, global_trajectory(:,1));
                %resample global_trajectory
                global_trajectory = interp1(selected_time, global_trajectory(:,1:3), selected_time(1):0.05:selected_time(end));

                %Resampling trajectory. 
                selected_time = selected_time(1):0.01:selected_time(end);

                f1 = figure;
                f1.Name = sprintf('%s - Animation',args.label);
                ax3 = axes(f1);hold on;axis(ax3,'equal');

                global_trajectory = global_trajectory(:,1:3);

                plot(path(:,1),path(:,2),'Parent',ax3,'LineStyle','--');
                
                x_range = [min(global_trajectory(:,1)),max(global_trajectory(:,1))];
                y_range = [min(global_trajectory(:,2)),max(global_trajectory(:,2))];
                x_lim = mean(x_range) + [-.5 .5]*max( [(x_range(2)-x_range(1)), (y_range(2)-y_range(1))]);
                y_lim = mean(y_range) + [-.5 .5]*max( [(x_range(2)-x_range(1)), (y_range(2)-y_range(1))]);
                y_lim(2) = y_lim(2) + 1e-10 * (y_lim(2) == y_lim(1));
                x_lim(2) = x_lim(2) + 1e-10 * (x_lim(2) == x_lim(1));
                set(ax3,'XLim',x_lim,'YLim',y_lim);

                plot(global_trajectory(:,1),global_trajectory(:,2),'Parent',ax3,'LineWidth',1);
                plot_handle = plot(global_trajectory(1,1),global_trajectory(1,2),'Parent',ax3,'LineWidth',2);
                
                animation_step_size = 0.03;
                n_steps = ceil(animation_step_size/mean(diff(selected_time)));
                car_profile = obj.carProf(parameters(2),0.5*parameters(2));
                car_profile(:,1) = car_profile(:,1) + 0.5*parameters(2);
                
                frame_cnt = 0;
                for i = 1:n_steps:(size(global_trajectory,1)-1)
                    rotated_car_profile = obj.rotatePath(car_profile,global_trajectory(i,3));
                    set(plot_handle,'xdata',rotated_car_profile(:,1) + global_trajectory(i,1),'ydata',rotated_car_profile(:,2) + global_trajectory(i,2));
                    
                    x_lim = global_trajectory(i,1) + [-1 1]*parameters(2)*5;
                    y_lim = global_trajectory(i,2) + [-1 1]*parameters(2)*5;
                    
                    if x_lim(2)<=x_lim(1), x_lim(2) = x_lim(1) + 0.0001;end
                    if y_lim(2)<=y_lim(1), y_lim(2) = y_lim(1) + 0.0001;end
                    if ~isnan(x_lim(1)) && ~isnan(y_lim(1)), set(ax3,'XLim',x_lim,'YLim',y_lim); end
                        
                    %saveImg(f1,sprintf('frame%05d.jpg',frame_cnt),'/Users/Ezequiel/Desktop/video_example/frames/',1,2)
                    frame_cnt = frame_cnt + 1;

                    drawnow;

                end
            end
        end
        %
        function path = rotatePath(~, path,ang)
            % path = rotatePath(PATH,ANG) takes a two columns array %
            % representing a path and rotates its by a certain angle ANG.
            rot = [cos(ang),sin(ang);-sin(ang),cos(ang)];
            path = path*rot;
        end
        %
        function compareObservedVsRealStates(obj, varargin)
            % compareObservedVsRealStates(TIME, DISCRETE_T, OBSERVED_X, NOISY_X, X) 
            % plots the real states of the system X, the noisy states X,
            % and the estimation of the states given by the implemented
            % observer OBSERVED_X in the same plot for comparison purposes.
            % The function requires the TIME and DISCRETE_TIME vector as
            % additional inputs.
            %
            % Additionally, those inputs can be followed by some
            % parameter/value pairs to specify additional properties of the
            % script.
            %
            %   - figure_name : string to be set as figure title.            
            p = inputParser;
            p.addRequired('time')
            p.addRequired('discrete_time')
            p.addRequired('observed_state')
            p.addRequired('noisy_state')
            p.addRequired('state')
            
            p.addOptional('figure_name',[])
            p.parse(varargin{:})
            args = p.Results;            
            
            ax_cells = cell(3,1);
            fig_handle = figure;
            n_states = size(args.state,2);
            ax_cells{2} = gobjects(n_states,1);
            
            label_text = 'State';
            for i = 1:n_states
                ax_cells{2}(i) = subplot(n_states,1,i,'Parent',fig_handle,'NextPlot','add');
                title(sprintf('%s - %d',label_text,i));
            end
            if ~isempty(args.figure_name), fig_handle.Name = sprintf('%s - states',args.figure_name); end
            obj.drawColumnsSeparately(args.time, args.discrete_time, args.noisy_state, ax_cells{2},'noisy state', 'x_label', 'time', 'y_labels', {'s [m]','d [m]','\theta_e [rad]','v [m/s]','\phi [rad]'}, 'plot_args',{'Color',[1 1 1]*.8,'LineWidth',0.5});
            obj.drawColumnsSeparately(args.time, args.discrete_time, args.observed_state, ax_cells{2}, 'observed state', 'plot_args', {'Color',[0 0 1],'LineWidth',2});
            obj.drawColumnsSeparately(args.time, args.discrete_time, args.state, ax_cells{2}, 'real state', 'plot_args', {'Color',[0 0 0],'LineWidth',1});
        end
        %
        function drawColumnsSeparately(obj, time, discrete_time, y_data, ax_cells, plot_label, varargin)
            % drawColumnsSeparately(TIME, DISCRETE_TIME, Y_DATA, AX_CELLS, PLOT_LABEL)
            % plots the columns of Y_DATA in the axes handle found within
            % the cell AX_CELLS. As Y_DATA is expected to be a simulation
            % output, the TIME and DISCRETE_TIME vectors are required.
            % Moreover, for the sake of clarity, a PLOT_LABEL must be
            % specified which will be shown in the legen and would identify
            % each plot. 
            %
            % Additionally, those inputs can be followed by some
            % parameter/value pairs to specify additional properties of the
            % script.
            %
            %   - noisy_y_data : must be specified if a noisy version of
            %   Y_DATA is wanted to be represented in the background. 
            %   - x_label : string to be set as X_LABEL.
            %   - y_labels : string to be set as Y_LABEL.
            %   - plot_args : cell which would include additional inputs to
            %   be pass to the plot function. 
            %   - stairs : boolean flag which would make this function use
            %   the stairs instead of the plot function to represent the
            %   given data.
            %
            p = inputParser;
            p.addOptional('noisy_y_data',[])
            p.addOptional('x_label',[])
            p.addOptional('y_labels',[])
            p.addOptional('plot_args',{})
            p.addOptional('stairs',true)
            p.parse(varargin{:})
            args = p.Results;            
            selected_time = obj.select_time_vector(time, discrete_time, y_data);
            
            for i = 1:size(y_data,2)
                if ~isempty(args.noisy_y_data)
                    noisy_t = obj.select_time_vector(time, discrete_time, args.noisy_y_data);
                    plot(noisy_t,args.noisy_y_data(:,i),'Parent',ax_cells(i),'color',[1 1 1]*.8);
                    line_width = 2;
                else
                    line_width =0.5;
                end
                    
                if args.stairs, stairs(selected_time, y_data(:,i),'Parent',ax_cells(i),'LineWidth',line_width, args.plot_args{:});
                else, plot(selected_time, y_data(:,i),'Parent',ax_cells(i),'LineWidth',line_width, args.plot_args{:});
                end
                if ~isempty(plot_label)
                    legend_texts = legend(ax_cells(i));
                    if ~isempty(legend_texts), legen_strings = legend_texts.String(1:(end));
                    else, legen_strings = {}; end
                    strings = {legen_strings{:}, plot_label};
                    valid_idx = ones(length(strings),1);
                    for k = 1:length(valid_idx)
                        if length(strings{k})>4 && strcmp(strings{k}(1:4),'data')
                            valid_idx(k) = 0;
                        end
                    end
                    strings = strings(valid_idx>0);
                    legend(ax_cells(i), strings{:});
                end
                if ~isempty(args.x_label), xlabel(ax_cells(i),args.x_label); end
                if ~isempty(args.y_labels), ylabel(ax_cells(i),args.y_labels{i}); end
            end
        end
        %
        function drawSingularValEvol(~, singular_values)
            % drawSingularValEvol(SINGULAR_VALUES) plots the evolution of the SINGULAR_VALUES to
            %represent the convergency of matrix S as the time horizon
            %increases.
            f = figure('Name','Evolution of svds of S for convergency illustration');
            ax = axes('Parent',f,'NextPlot','add');
            plot(ax, singular_values);
            grid;
            xlabel('Nombre of iterations'),ylabel('singular velues of matrix S')            
        end
        %
        function [Nodes, Edges] = buildDynProgGraph(obj, sampling_time, initial_state, path_to_track, allowed_u, parameters, varargin)
            % [NODES, EDGES] = buildDynProgGraph(SAMPLING_TIME, INITIAL_STATE,
            % PATH_TO_TRACK, ALLOWED_U, PARAMETERS) creates a graph
            % (defined by a set of NODES and EDGES) representing the
            % possible trajectories that the system can follow given its
            % INITIAL_STATE, the set of control inputs ALLOWED_U, the
            % PARAMETERS of the system, the PATH_TO_TRACK and the %
            % SAMPLING_TIME between control actions. 
            %
            % In order to alleviate the computional complexity required to
            % build a fully connected tree, a randomized approach is herein
            % adopted. Specifically, once a set of fully connected nodes
            % are built, the number of new nodes that are created from %
            % every parent node, will be dynamically reduced. 
            %
            %
            % This function returns two outputs: 
            %  Nodes = [ID, t, flag, px, py, yaw, x1 x2 x3 x4 x5] which is
            %  an array with as many rows as nodes in the graph, containing
            %  information regarding the state that each node represents as
            %  well as the TIME, NODEID and an additional FLAG. 
            %  - ID: mode ID
            %  - t:  time associated to the node
            %  - flag: if the flag is:
            %       - 0 if the node has not been visited. This is the case
            %       of the last node in the branches of the tree. 
            %       - 1 if the node was visited. 
            %       - 2 if the node is a dead end. 
            %  - [px, py, yaw] : triple showing the absolute position of
            %  the vehicle. 
            %  - [x1 x2 x3 x4 x5] : state of the vehicle that node
            %  represents.
            %
            %
            % The ouput array EDGES contains the information of which paris 
            % of nodes are connected. Specifically it will be a 4-columns
            % array
            % containing, in order, 
            %  - predecessor node id
            %  - sucessor node id. 
            %  - [u1,u2] control input applied between nodes. 
            %
            % Additionally, those inputs can be followed by some
            % parameter/value pairs to specify additional properties of the
            % script.
            %
            %   - max_n_nodes : maximum number of nodes to be included in
            %   the graph.  
            %   - plot_nodes : boolean flag showing whether the nodes are
            %   wanted to be represented.
            %   - plot_edges : boolean flag showing wether the nodes are
            %   wanted to be represented. 
            %   - plot : boolean flag showing whether the full graph should
            %   be represented.
            %   - max_heading_deviation : max heading deviation a node must
            %   have so that it is considered innadmisible and therefore
            %   flaged as dead end. 
            %   - max_lateral_deviation : max lateral deviation a node must
            %   have so that it is considered innadmisible and therefore
            %   flaged as dead end. 
            %   - time_to_loose_it : time indicating the time from which
            %   the branches of the tree will be randomly reduced. 
            %   - threshold_number_of_nodes : minimum number of nodes that
            %   are wanted to be included for every time. The implemented
            %   randomized branch construction will use this parameter to
            %   adapt the number of edges to be constructed from every node
            %   in the previous layer. 
            p = inputParser;
            p.addOptional('max_n_nodes',25000);
            p.addOptional('plot_nodes',false);
            p.addOptional('plot_edges',false);
            p.addOptional('plot',false);
            p.addOptional('max_heading_deviation',pi/4);
            p.addOptional('max_lateral_deviation',3.5);
            p.addOptional('time_to_loose_it',2);
            p.addOptional('threshold_number_of_nodes',120);
            
            p.parse(varargin{:});
            args = p.Results;
            
            n_states = length(initial_state) + 3;
            % initializing arrays
            Nodes = zeros(args.max_n_nodes, 3 + n_states);
            Edges = zeros(args.max_n_nodes, 4);
            
            %Setting the nodes ID. 
            Nodes(:,1) = (1:args.max_n_nodes)';
            
            %set up first node. 
            path_ori =  interp1(path_to_track.s_coordinate, path_to_track.orientation, initial_state(1));
            initial_position = ...
                interp1(path_to_track.s_coordinate, path_to_track.coordinates, initial_state(1)) + ...
                [-sin(path_ori), cos(path_ori)]*initial_state(2);
            
            orientation = path_ori + initial_state(3); 
            Nodes(1,4:end) = [initial_position,orientation,initial_state'];
            
            next_node_pointer = 2;
            next_edge_pointer = 1;
            end_flag = false;
            
            cnt = 0;
            tic
            fprintf('Building the graph');
            while ~end_flag 
                % get the list of nodes that should be pushed forward. 
                open_nodes_id_list = find( Nodes(:,3) == 0 & Nodes(:,1) < next_node_pointer);
                
                % for every of those nodes
                for root_node_id = open_nodes_id_list'
                    cnt = obj.printProgression(cnt, args.max_n_nodes-args.threshold_number_of_nodes, 100);
                    valid = false;
                    
                    % for every u option
                    u_options = 1:size(allowed_u,1);
                    while ~valid && any(u_options>0) && ~end_flag
                        
                        chosen_allowed_u = allowed_u;
                        chosen_idxs = 1:size(chosen_allowed_u,1);
                        valid_idxs = 1:size(chosen_allowed_u,1);
                        
                        %retrieve the state of the root node. 
                        root_node_state = Nodes(root_node_id,4:end);
                        
                        %how many nodes in the previous instant were dead?
                        nodes_ID_prev_time = find((Nodes(:, 2) == Nodes(root_node_id, 2)));
                        n_nodes_previous_time = length(nodes_ID_prev_time);
                        previous_time_dead_nodes = sum( (Nodes(:, 2) == Nodes(root_node_id, 2)) & (Nodes(:, 3) == 2));
                        
                        randomize = Nodes(root_node_id, 2)>args.time_to_loose_it;
                        
                        enough_nodes_prev_time = sum( abs(Nodes(nodes_ID_prev_time,8))<2) >args.threshold_number_of_nodes;
                        enough_nodes_prev_time = (n_nodes_previous_time - previous_time_dead_nodes)>args.threshold_number_of_nodes;
                        
                        if randomize && enough_nodes_prev_time
                            %chose only u_options control values to
                            %evaluate. 
                            
                            u_ops = find(u_options);
                            if root_node_state(6)>0
                                u_ops = u_ops( allowed_u(u_ops,2)<=0);
                            else
                                u_ops = u_ops( allowed_u(u_ops,2)>=0);
                            end
                            if isempty(u_ops)
                                u_ops = find(u_options);
                            end
                            
                            
                            valid_idxs = find(u_ops>0);
                            valid_idxs = u_ops;
                            valid_idxs = find(u_options);
                            chosen_idxs = randi(length(valid_idxs), min(length(valid_idxs), 1),1);
                                                        
                            chosen_allowed_u = allowed_u(valid_idxs(chosen_idxs),:);
                        end
                        
                        %id of the following nodes, which we now because we now
                        %how many actions combinations we are evaluating
                        n_new_states = size(chosen_allowed_u,1);
                        new_nodes_id = next_node_pointer + (1:n_new_states) -1;
                        chosen_allowed_u = chosen_allowed_u(new_nodes_id<args.max_n_nodes,:);
                        new_nodes_id = new_nodes_id(new_nodes_id<args.max_n_nodes);
                        
                        n_new_states = length(new_nodes_id);




                        root_node_states = repmat(root_node_state',1,n_new_states);
                        new_states = root_node_states;

                        %curvature = interp1(path_to_track.s_coordinate, path_to_track.curvature, rem(root_node_states(4,:),path_to_track.s_coordinate(end)));

                        %non linear model implementation
                        new_states(1,:) = root_node_states(1,:) + sampling_time*( root_node_states(7,:) .* cos(root_node_states(3,:)) );
                        new_states(2,:) = root_node_states(2,:) + sampling_time*( root_node_states(7,:) .* sin(root_node_states(3,:)) );
                        new_states(3,:) = root_node_states(3,:) + sampling_time*( root_node_states(7,:)/parameters(2).* tan(root_node_states(8,:)) );

                        dist_path_to_states = permute(sum((path_to_track.coordinates - permute(new_states(1:2,:)',[3,2,1])).^2,2),[1,3,2]);
                        min_dist_to_path = min(dist_path_to_states,[],1);
                        proj_point_id = sum((dist_path_to_states == min_dist_to_path) .* ((1:size(dist_path_to_states,1))'),1);
                        s_proj = path_to_track.s_coordinate(proj_point_id);
                        head = path_to_track.orientation(proj_point_id);
                        
                        %ds =  (root_node_states(7,:).*cos(root_node_states(6,:)))./(1 - root_node_states(6,:).*curvature);
                        %new_states(4,:) = root_node_states(4,:) + sampling_time*( ds );
                        %new_states(5,:) = root_node_states(5,:) + sampling_time*( root_node_states(7,:).*sin(root_node_states(6,:)) );
                        %new_states(6,:) = root_node_states(6,:) + sampling_time*( root_node_states(7,:)./parameters(2).*tan(root_node_states(8,:)) - curvature .* ds );

                        
                        new_states(4,:) = s_proj;
                        new_states(5,:) = min(dist_path_to_states,[],1);
                        new_states(6,:) = new_states(3,:) - head';                        
                        
                        new_states(7,:) = root_node_states(7,:) + sampling_time*(parameters(3) * (chosen_allowed_u(:,1)' - root_node_states(7,:)));
                        new_states(8,:) = root_node_states(8,:) + sampling_time*(parameters(4) * (chosen_allowed_u(:,2)' - root_node_states(8,:)));                  


                        if any(new_nodes_id)>args.max_n_nodes
                            error('You are changing the size of Nodes array, there is something wrong');
                        end
                        Nodes(new_nodes_id,2) = Nodes(root_node_id,2) + sampling_time;

                        Nodes(new_nodes_id,4:end) = new_states';

                        %set the root node as explored
                        Nodes(root_node_id,3) = 1;

                        %build the edges
                        Edges( next_edge_pointer + (1:n_new_states)-1,1) = root_node_id;
                        Edges( next_edge_pointer + (1:n_new_states)-1,2) = new_nodes_id;
                        Edges( next_edge_pointer + (1:n_new_states)-1,3:4) = chosen_allowed_u;

                        next_edge_pointer = next_edge_pointer + n_new_states;

                        %if some of the new nodes are already too far away,
                        %mark them as explored even if they are not, so that
                        %they do not keep on being broadcasted. 
                        node_state = Nodes(new_nodes_id, 7:11);
                        Nodes(new_nodes_id, 3) = 2*( abs(node_state(:, 2)) > args.max_lateral_deviation );
                        Nodes(new_nodes_id, 3) = 2*( Nodes(new_nodes_id, 3)==2 | (abs(node_state(:, 3)) > args.max_heading_deviation));
                        Nodes(new_nodes_id, 3) = 2*( Nodes(new_nodes_id, 3)==2 | (abs(node_state(:, 4)) < 0));

                        next_node_pointer = new_nodes_id(end) + 1;
                        
                        valid = any(Nodes(new_nodes_id, 3)==0);
                        u_options(valid_idxs(chosen_idxs)) = -1*(Nodes(new_nodes_id, 3)==1);
                        %if next_node_pointer > args.max_n_nodes
                        %    ''
                        %end
                        end_flag = end_flag | (next_node_pointer > (args.max_n_nodes - 3));
                    end
                end
                end_flag = end_flag | isempty(open_nodes_id_list);
            end
            fprintf('\n');
            if args.plot_nodes || args.plot_edges || args.plot
                figure;
                ax = axes('NextPlot','add');
                n_valid_edges = find(Edges(:,1)==0,1,'first') - 1;
                
                plot(path_to_track.coordinates(:,1),path_to_track.coordinates(:,2),'Parent',ax)
                if args.plot_nodes && ~args.plot_edges && ~args.plot
                    xdata = Nodes(:,4);
                    ydata = Nodes(:,5);                    
                    plot(xdata, ydata,'.','MarkerSize',5,'Parent',ax,'Color',[0 0 0]);
                elseif ~args.plot_nodes && args.plot_edges && ~args.plot
                    xdata = [Nodes(Edges(1:n_valid_edges,1),4), Nodes(Edges(1:n_valid_edges,2),4)]';
                    ydata = [Nodes(Edges(1:n_valid_edges,1),5), Nodes(Edges(1:n_valid_edges,2),5)]';
                    plot(xdata,ydata,'MarkerSize',5,'Parent',ax,'Color',[0 0 0]);
                else
                    xdata = [Nodes(Edges(1:n_valid_edges,1),4), Nodes(Edges(1:n_valid_edges,2),4)]';
                    ydata = [Nodes(Edges(1:n_valid_edges,1),5), Nodes(Edges(1:n_valid_edges,2),5)]';
                    plot(xdata,ydata,'.-','MarkerSize',5,'Parent',ax,'Color',[0 0 0]);
                end
                
                
            end
            
            
            Nodes = Nodes( ~ (Nodes(:,3) == 0 & Nodes(:,2) == 0) ,:);
            Edges = Edges(Edges(:,1)~=0,:);
        end
        %
        function counter = printProgression(~,  counter, n_iterations, n_slots)
            % printProgression(COUNTER, N_ITERATIONS, N_SLOTS) prints a
            % loading bar in the command windows. It takes the number of
            % iterations N_ITERATIONS that are expected, the current
            % iteration COUNTER and the number of slots N_SLOTS that are
            % wanted to be represented. It then returns the updated state
            % of the iteration COUNTER.
            persistent n_printed_intervals
            if isempty(n_printed_intervals) 
                n_printed_intervals = 0;
            end
            n_slots = min(n_slots, n_iterations-2);
            if counter == 0
                n_printed_intervals = 0;
                fprintf(['\n|' repmat(' ',1,n_slots-2) '|\n ']);
            end
            iterations_interval = ceil(n_iterations/n_slots);
            %n_intervals_already_printed = max(0,floor((cnt-1)/iterations_interval));
            %if cnt == 0 || (cnt-iterations_interval*floor(cnt/iterations_interval-0.00000001) )>=floor(iterations_interval)
            if counter == 0 || (counter-iterations_interval*n_printed_intervals )>=iterations_interval
                n_printed_intervals = n_printed_intervals + 1;
                fprintf('*');
            end
            if counter>=(n_iterations-1), fprintf('*\n'); end
            counter = counter + 1;
        end
        %
        function [ax] = plotGraph(~, nodes, edges, path_to_track, varargin)
            % plotGraph(NODES, EDGES, PATH_TO_TRACK) draws the graph 
            % branches and ell as the PATH_TO_TRACK for comparison purposes. 
            %
            % Additionally, those inputs can be followed by some
            % parameter/value pairs to specify additional properties of the
            % script.
            %
            %   - PLOT_NODES : boolean flag that shows whether the nodes
            %   have to be represented. 
            %   - PLOT_EDGES : boolean flag that shows whether the edges
            %   have to be represented. 
            %   - PLOT : boolean flag showing whether the entire graph
            %   should be represented.
            %   - PLOT_TERMINAL_NODES : boolean flagh showing whether the
            %   terminal nodes should be highlighted. 
            %   - PLOT_DEAD_ENDS : boolean flag showing whether the nodes
            %   representing dead ends should be highlighted.
            %   - obstacles : structure with a unique field 'position'
            %   containing a two-columns array (with as many rows as 
            %   obstacles are wanted to be included) representing the %
            %   position of a set of obstacles.
            %   - obstacle_marker_size : parameter showing the size of the
            %   marker representing the obstacles. 
            p = inputParser;            
            p.addOptional('plot_nodes',false);
            p.addOptional('plot_edges',false);
            p.addOptional('plot',false);
            p.addOptional('plot_terminal_nodes',false);
            p.addOptional('plot_dead_ends',false);
            p.addOptional('obstacles',[]);
            p.addOptional('obstacle_marker_size',40);
            
            p.parse(varargin{:});
            args = p.Results;
            
            figure;
            ax = axes('NextPlot','add');
            n_valid_edges = find(edges(:,1)==0,1,'first') - 1;
                
            plot(path_to_track.coordinates(:,1),path_to_track.coordinates(:,2),'Parent',ax)            
            
            if args.plot_nodes || args.plot_edges || args.plotargs.plot
                nodes_color = [1 1 1]*.7;
                if args.plot_nodes && ~args.plot_edges && ~args.plot
                    xdata = nodes(:,4);
                    ydata = nodes(:,5);                    
                    plot(xdata, ydata,'.','MarkerSize',5,'Parent',ax,'Color',nodes_color);
                elseif ~args.plot_nodes && args.plot_edges && ~args.plot
                    xdata = [nodes(edges(1:n_valid_edges,1),4), nodes(edges(1:n_valid_edges,2),4)]';
                    ydata = [nodes(edges(1:n_valid_edges,1),5), nodes(edges(1:n_valid_edges,2),5)]';
                    plot(xdata,ydata,'MarkerSize',5,'Parent',ax,'Color',nodes_color);
                else
                    xdata = [nodes(edges(1:n_valid_edges,1),4), nodes(edges(1:n_valid_edges,2),4)]';
                    ydata = [nodes(edges(1:n_valid_edges,1),5), nodes(edges(1:n_valid_edges,2),5)]';
                    plot(xdata,ydata,'.-','MarkerSize',5,'Parent',ax,'Color',nodes_color);
                end
            end
            
            if args.plot_terminal_nodes
                terminal_nodes = find(nodes(:,3)==0);
                xdata = nodes(terminal_nodes,4);
                ydata = nodes(terminal_nodes,5);
                plot(xdata, ydata,'.','MarkerSize',20,'Parent',ax,'Color',[0 0 1]);                
            end
            
            if args.plot_dead_ends
                dead_end_nodes = find(nodes(:,3)==2);
                xdata = nodes(dead_end_nodes,4);
                ydata = nodes(dead_end_nodes,5);
                plot(xdata, ydata,'.','MarkerSize',10,'Parent',ax,'Color',[1 0 0]);
            end          
            
            if ~isempty(args.obstacles) && ~isempty(args.obstacles.position)
                plot(args.obstacles.position(:,1), args.obstacles.position(:,2),'.','MarkerSize',args.obstacle_marker_size,'Parent',ax,'Color',[0.5 0 0]);
            end
        end
        %
        function plotDynProgTraj(~,  nodes, optimal_node_sequence, ax)
            % plotDynProgTraj(NODES, OPTIMAL_NODE_SEQUENCE, AX)
            % represents in the axes AX, the path defined by the
            % OPTIMAL_NDOE_SEQUENCE, by using the information in NODES. 
            if ~isempty(ax)
                xdata = nodes(optimal_node_sequence(1:(end-1)), 4);
                ydata = nodes(optimal_node_sequence(1:(end-1)),5);
                plot(xdata,ydata,'Color',[0 0 0],'Parent',ax,'LineWidth',2);
            end
        end
        %
        function car = carProf(~, car_length, car_width)
            % CAR = carProf(CAR_LENGTH, CAR_WIDTH) returns a two-columns
            % array defining the profile of a car given its CAR_LENGTH and
            % CAR_WIDTH, which will be used for representation purposes. 
            car =  [-0.4537,-0.4081,-0.3727,-0.3473,-0.3320,-0.3224,-0.3143   -0.3078   -0.3029   -0.2995   -0.2934   -0.2823   -0.2688   -0.2555   -0.2453   -0.2402   -0.2395   -0.2421   -0.2466   -0.2518   -0.2537   -0.2502   -0.2431   -0.2344   -0.2258   -0.1816   -0.0867    0.0249    0.1197    0.1640    0.1736    0.1866 0.2028    0.2221    0.2441    0.2694    0.2968    0.3237    0.3474    0.3654    0.3900    0.4297    0.4761    0.5210    0.5560    0.5958    0.6536    0.7188    0.7808    0.8289    0.8708    0.9163    0.9579    0.9882    1.0000    0.9907    0.9581    0.8947    0.7935    0.6470    0.5183    0.4548    0.4252    0.3986 0.3438    0.2841    0.2544    0.2412    0.2307    0.2095    0.1864    0.1731    0.1604    0.1388    0.0990    0.0320   -0.0538   -0.1411   -0.2128   -0.2518   -0.2655   -0.2731   -0.2754   -0.2731   -0.2670   -0.2589   -0.2520   -0.2492   -0.2533   -0.2670   -0.2846   -0.2993   -0.3124   -0.3250   -0.3385   -0.3594 -0.3909   -0.4291   -0.4699   -0.5096   -0.5711   -0.6633   -0.7606   -0.8377   -0.8691   -0.8415   -0.7719   -0.6801   -0.5860   -0.5096   -0.4537;
                    0.7740    0.7680    0.7601    0.7520    0.7453    0.7430    0.7460    0.7533    0.7638    0.7764    0.7955    0.8219    0.8490    0.8700    0.8785    0.8698    0.8467    0.8143    0.7774    0.7409    0.7159    0.7065    0.7070    0.7114    0.7142    0.7120    0.7073    0.7027    0.7011    0.7053    0.7154    0.7269 0.7360    0.7390    0.7320    0.7202    0.7114    0.7054    0.7020    0.7009    0.7046    0.7134    0.7239    0.7327    0.7364    0.7373    0.7374    0.7326    0.7188    0.6920    0.6518    0.5821    0.4557    0.2454   -0.0763   -0.4136   -0.6593   -0.8229   -0.9140   -0.9423   -0.9353   -0.9188   -0.8991   -0.8826 -0.8757   -0.8816   -0.8957   -0.9121   -0.9251   -0.9289   -0.9222   -0.9101   -0.8970   -0.8870   -0.8845   -0.8834   -0.8762   -0.8657   -0.8543   -0.8446   -0.8394   -0.8404   -0.8484   -0.8643   -0.8890   -0.9206   -0.9530   -0.9807   -0.9983   -1.0000   -0.9848   -0.9601   -0.9338   -0.9135   -0.9067   -0.9103 -0.9149   -0.9195   -0.9231   -0.9245   -0.9271   -0.8986   -0.7845   -0.5300   -0.0807    0.3700    0.6277    0.7458    0.7776    0.7764    0.7740]';
            car(:,1) = -car(:,1) * 0.5*car_length;
            car(:,2) = car(:,2) * 0.5*car_width;
        end
        %
        function sol_class = load_solution_class(obj,  exercise_id)
            % SOL_CLASS = load_solution_class(EXERCISE_ID) allows to load
            % the solution class for exercise EXERCISE_ID.
            master_computer = obj.isMasterComputer;
            switch exercise_id
                case 1, if master_computer, sol_class = solEx1; else, sol_class = ex1; end
                case 2, if master_computer, sol_class = solEx2; else, sol_class = ex2; end
                case 3, if master_computer, sol_class = solEx3; else, sol_class = ex3; end
                case 4, if master_computer, sol_class = solEx4; else, sol_class = ex4; end
            end
        end
        %
        function complete_name = complete_simulink_model_name(obj,  name_simulink_file)
            % COMPLETE_NAME =
            % complete_simulink_model_name(name_simulink_file) returns a
            % string containing the path of the right simulink file to use
            % depending on your matlab version and the NAME_SIMULINK_FILE.
            % It requires the provided directory structure to stay unchanged. 
            %matlab_version = version('-release');
            extension = '.slx';
            if obj.isMasterComputer
                solved_label = 'solved_';
            else
                solved_label = '';
            end
            if isfolder('./functions/super_user_utilities/')
                complete_name = [solved_label, name_simulink_file, extension];
            else
                file = dir(['./', solved_label, name_simulink_file, '*', extension]);
                complete_name = file.name;
            end
        end
        %
        function flag = isMasterComputer(~)
            % FLAG = isMasterComputer return a booelan flag if the
            % computer where the script is being run is the teacher's
            % computer. This function changes the behavior of a few
            % functions within this utility class, mainly related to the
            % solution classes and simulink files that are loaded. Notice
            % that changing manually this value will not magically provide
            % you with the solution of the exercise but rather braking the 
            % code.  
            flag = isfile('key.txt');
        end
        %
        function printLabeledAB(obj,  varargin)
            % printLabeledAB(MATRIX_A, MATRIX_B, LABEL) prints in the
            % command window the value of matrices MATRIX_A and MATRIX_B in
            % a clean way and preceded by the string pass within LABEL.    
            p = inputParser;
            p.addRequired('matrix_A');
            p.addRequired('matrix_B');
            p.addRequired('label');
            
            p.addOptional('description', '[Phi | Gamma]');
            p.parse(varargin{:});
            args = p.Results;
            
            fprintf('\n%s = %s =\n', args.label, args.description);
            obj.printJointAB( args.matrix_A, args.matrix_B);
        end
        %
        function printRelativeError(obj, varargin)
            % printRelativeError(MATRIX_A1, MATRIX_A2, MATRIX_B1,
            % MATRIX_B2, MATRIX_LABEL) prints in the command window, in a
            % clean way, the relative error between MATRIX_A1, MATRIX_A2,
            % and MATRIX_B1, MATRIX_B2. The numerical values of the
            % relative error is preceded by the string in MATRIX_LABEL
            p = inputParser;
            p.addRequired('matrix_A1');
            p.addRequired('matrix_A2');
            p.addRequired('matrix_B1');
            p.addRequired('matrix_B2');
            p.addRequired('matrix_label');
            
            p.addOptional('description', '[relative error Phi | relative error Gamma ]');
            p.parse(varargin{:});
            args = p.Results;
            
            errorA = abs(args.matrix_A1-args.matrix_A2)/(1e-10 + args.matrix_A1);
            errorB = abs(args.matrix_B1-args.matrix_B2)/(1e-10 + args.matrix_B1);
            obj.printLabeledAB( errorA, errorB, args.matrix_label, 'description', args.description);
        end        
        %
        function printJointAB(~, Phi,Gam)
            % printJointAB(PHI, GAM) prints in the command window the PHI
            % and GAM matrices jointly. 
            for i = 1:size(Phi,1)
                fprintf('\t[');
                fprintf('%8.3f ', Phi(i,:));
                fprintf(' | ');
                fprintf('%8.3f ', Gam(i,:));
                fprintf(']\n');
            end            
        end
        %
        function noiseMod = defaultNoiseModule(~)
            % noiseMod = defaultNoiseModule returns the noise structure
            % used to introduce noise in the simulations with some default
            % values. 
            noiseMod.mean = [0; 0; 0; 0; 0];
            noiseMod.std_deviation = [0; 0; 0; 0; 0];
            noiseMod.seed = [0; 0; 0; 0; 0];
            noiseMod.freq = 5;
        end
        %
        function printComplexNumbers(~, data, varargin)
            % printComplexNumbers(DATA,LABEL) prints in the command windows and in a
            % clean way, the set of complex number passed in DATA.
            % Information that should be preceded by the string in LABEL. 
            p = inputParser;
            p.addOptional('label','');
            p.parse(varargin{:});
            args = p.Results;
            
            fprintf('\n%s : \n', args.label)
            R = real(data);
            IM = imag(data);
            for j = 1:length(data)
                fprintf('\t%10.7f, %10.7f i\n',R(j),IM(j))
            end
        end
        %
        function stop_flag = stopCondition(~,  function_output, function_name)
            % SIMULATE_FLAG = stopCondition(FUNCTION_OUTPUT,
            % FUNCTION_NAME) returns a boolean variable
            % SIMULATE_FLAG which is false if the FUNCTION_OUTPUT given 
            % as input is empty.
            %
            % Additionally, this function prints a message pointing out the
            % functions that have not been completed, for which the input
            % string FUNCTION_NAME is used.
            stop_flag = false;
            if isempty(function_output)
                fprintf('\nComplete function %s to continue\n', function_name);
                stop_flag = true;
            end
        end
        %
        function XY_position = SDCoordToXY(~,  path, s_coord, d_coord)
            % XY_POSITION = SDCoordToXY( PATH_TO_TRACK, S_COORD, D_COORD)
            % returns the XY_POSITION of a path defined by the S_COORD and
            % D_COORD coordinates, and the PATH w.r.t. which they are
            % calculated.
            
            XY_position = interp1( path.s_coordinate, path.coordinates, s_coord);
            orientation = interp1( path.s_coordinate, path.orientation, s_coord);
            
            if size(s_coord,2)>1
                XY_position = XY_position + repmat(d_coord,1,1,2).*cat(3,cos(orientation + pi/2), sin(orientation + pi/2));
            else
                XY_position = XY_position + repmat(d_coord,1,2).*[cos(orientation + pi/2), sin(orientation + pi/2)];
            end
        end
        %
        function [ax, potential_field, s, d, s_grid, d_grid] = representPotentialFields(ut, path_to_track, attractive_field, repulsive_field, varargin)
            % [AX, POTENTIAL_FIELD, S, D, S_GRID, D_GRID] =
            % representPotentialFields(UT, PATH_TO_TRACK, ATTRACTIVE_FIELD,
            % REPULSIVE_FIELD) represents the artificial potential field
            % resulting from combining the ATTRACTIVE_FIELD and
            % REPULSIVE_FIELD passed as input, as well as the
            % PAHT_TO_TRACK. 
            p = inputParser;
            p.addOptional('plot',true);
            p.addOptional('s_step',1);
            p.addOptional('d_bounds',15);
            p.addOptional('d_step',0.2);
            p.parse(varargin{:});
            args = p.Results;    
                        
            s = 0:args.s_step:max(path_to_track.s_coordinate);
            d = -args.d_bounds:args.d_step:args.d_bounds;
            [s_grid, d_grid] = meshgrid(s,d);
            
            pos_grid = ut.SDCoordToXY( path_to_track, s_grid, d_grid);
            potential_field = -attractive_field(s_grid, d_grid) + repulsive_field(s_grid, d_grid);
            
            ax = [];
            if args.plot
                ax = gobjects(8,1);
                figure, ax(1) = axes; surf(attractive_field(s_grid, d_grid));
                figure;ax(2) = axes;surf(repulsive_field(s_grid, d_grid));             
                figure,ax(3) = axes; surf(attractive_field(s_grid, d_grid));
                figure,ax(4) = axes('NextPlot','add'); surf(s_grid, d_grid, potential_field);
                figure;ax(5) = axes('NextPlot','add');contour(pos_grid(:,:,1), pos_grid(:,:,2), potential_field,100,'Parent',ax(5))
                figure,ax(6) = axes;surf(pos_grid(:,:,1), pos_grid(:,:,2), attractive_field(s_grid, d_grid))
                figure,ax(7) = axes;surf(pos_grid(:,:,1), pos_grid(:,:,2), repulsive_field(s_grid, d_grid))
                figure,ax(8) = axes('NextPlot','add');surf(pos_grid(:,:,1), pos_grid(:,:,2), potential_field);
            end
        end
        %
        function path_s_d = obtainMinimumPotentialPath(~,initial_pos_s_d, s, d, fs, fd, varargin)
            % [PATH_S_D] = obtainMinimumPotentialPath(INITIAL_POS_S_D, S,
            % D, FS, FD) returns the PATH_S_D following the maximum
            % negative gradient of the potential field resulting from 
            % combining the input ATTRACTIVE_FIELD and REPULSIVE_FIELD. 
            % THe function also requires the S, D, vectors showing the grid 
            % used to generate the field terms as well as the gradient FS 
            % and FD of the potential field.           
            p = inputParser;
            p.addOptional('path_points',8000);
            p.addOptional('step_size',0.01);
            p.parse(varargin{:});
            args = p.Results;              

            gradient_data = zeros(args.path_points,2);
            path_s_d = zeros(args.path_points,2);
            path_s_d(1,:) = initial_pos_s_d;
            for counter = 2:args.path_points
                [~,  col] = min( abs(s - path_s_d(counter-1,1)));
                [~,  row] = min( abs(d - path_s_d(counter-1,2)));

                gradient_s_d = [ fs(row, col), 2 * fd(row, col) ];

                gradient_data(counter, :) = gradient_data(counter-1, :) + gradient_s_d;
                path_s_d(counter,:) = path_s_d(counter-1,:) - args.step_size*gradient_s_d;
            end
        end        
        %
        function circular_path = genCircularPath(~,  curvature)
            % circ_path = genCircularPath( curvature)
            % 
            % Returns a structure *circ_path* representing a circular path of
            % curvature *curvature* with fields:
            %   circ_path.coordinates   : Two-columns array [x, y];
            %   circ_path.s_coordinates : Column vector [s]
            %   circ_path.orientation   : Column vector [ori]
            %   circ_path.curvature     : Column vector [curvature]
            radius = 1/curvature;       
            theta_vector = (-pi/2+(0:0.0001:2*pi))';
            circular_path.coordinates = radius * [cos(theta_vector), sin(theta_vector)] + [zeros(size(theta_vector,1),1), ones(size(theta_vector,1),1)*radius];
            circular_path.s_coordinate = (theta_vector-theta_vector(1))*radius;
            circular_path.orientation = theta_vector+ pi/2;
            circular_path.curvature = circular_path.orientation *0 + curvature;            
        end
    end
    
end