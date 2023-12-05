%% Consensus with a static, undirected topology & F malicious agents
% MECH 6V29: MARS final project

% Demonstrates a working consensus algorithm where a given number of robots
% will experience failure & remain stationary for a undetermined amount of
% time
clear;clc;

N = 20;
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true);


%% Experiment constants

F = 3;                          % number of malicious agents
G = 3;                          % how many neighbors to ignore
iterations = 1000;              % number of iterations the experiment will run over
mal_r = randperm(N,F);          % list of which robots are malicious

% try out different topologies (make sure to change algorithm section
% accordingly)
% L = cycleGL(N);     % constant cycle graph
L = completeGL(N);  % constant complete graph
% L = ERGL(N,0.5);    % constant ER random graph


%% Grab tools we need to convert from single-integrator to unicycle dynamics

% Gain for the diffeomorphism transformation between single-integrator and
% unicycle dynamics
[~, uni_to_si_states] = create_si_to_uni_mapping();
si_to_uni_dyn = create_si_to_uni_dynamics_with_backwards_motion();

uni_barrier_cert_boundary = create_uni_barrier_certificate_with_boundary();

% Initialize velocity vector for agents.  Each agent expects a 2 x 1
% velocity vector containing the linear and angular velocity, respectively.
dxi = zeros(2, N);

%Create red marker for bad agents
marker_size_robot = determine_robot_marker_size(r);
CM = 'red';
line_width = 3;
font_size = determine_font_size(r, 0.05);

%plot iterations on bottom
iteration_caption = sprintf('Iteration %d', 0);
iteration_label = text(-1.5, -0.8, iteration_caption, 'FontSize', font_size, 'Color', 'r', 'FontWeight', 'bold');


%Iterate for the previously specified number of iterations
for t = 1:iterations

    % Retrieve the most recent poses from the Robotarium.  The time delay is
    % approximately 0.033 seconds
    x = r.get_poses();

    xi = uni_to_si_states(x);   % Convert to SI states
    xm = xi;                    % copy states to be corrupted

    % get corrupted states
    for state = 1:N
        for mal = 1:F
            if state == mal_r(mal)
                xm(:,state) = [1.6;1] - [3.2;2].*rand(2,1); % random (feasible) malicious value
            end
        end
    end

    %Initialize Bad Robot Markers
    for i = 1:N
        for robot = 1:F
            if i == mal_r(robot)
                % Plot colored circles showing robot location.
                g(i) = plot(x(1,i),x(2,i),'.','MarkerSize', marker_size_robot,'LineWidth',5,'Color',CM);
            end
        end
    end


    %% Algorithm

    for i = 1:N
        % Initialize velocity to zero for each agent.  This allows us to sum
        % over agent i's neighbors
        dxi(:, i) = [0 ; 0];

        malicious = false;  % initialize robot as non-malicious

        for robot = 1:F     % loop over total number of malicious nodes

            if i == mal_r(robot)    % robot is malicious
                malicious = true;   % indicate maliciousness
            end

        end


        % % do not update velocity if robot is malicious
        % if malicious
        %     continue
        % end

        % static graph
        neighbors = topological_neighbors(L, i);    % get list of neighbors
        distance = zeros(1,length(neighbors));      % initialize distance to neighbors

        % find (malicious) distance to neighbors
        for nhbr = 1:length(neighbors)
            distance(nhbr) = norm(xm(:,i) - xm(:,neighbors(nhbr)));
        end

        % ensure we have the confidence level to get rid of information
        % otherwise robot will not move
        if length(neighbors) > G

            [sorted_dist,nhbr_order] = sort(distance);  % sort distances

            % get rid of G furthest away neighbors & apply consnesus dynamics
            for j = 1:length(sorted_dist)-G
                dxi(:, i) = dxi(:, i) + (xm(:,neighbors(nhbr_order(j))) - xm(:, i));
            end

        end
    end
    % Update Iteration
    iteration_caption = sprintf('Iteration %d', t);
    iteration_label.String = iteration_caption;

    %% Avoid actuator errors
    
    % To avoid errors, we need to threshold dxi
    norms = arrayfun(@(x) norm(dxi(:, x)), 1:N);
    threshold = 3/4*r.max_linear_velocity;
    to_thresh = norms > threshold;
    dxi(:, to_thresh) = threshold*dxi(:, to_thresh)./norms(to_thresh);
    
    %% Map SI to Uni dynamics and utilize barrier certificates
    
    % Transform the single-integrator to unicycle dynamics using the the
    % transformation we created earlier
    dxu = si_to_uni_dyn(dxi, x);
    
    dxu = uni_barrier_cert_boundary(dxu, x);
    
    %% Send velocities to agents
     
    % Set velocities of agents 1,...,N
    r.set_velocities(1:N, dxu);
    
    % Send the previously set velocities to the agents.  This function must be called!    
    r.step();
end

% We can call this function to debug our experiment!  Fix all the errors
% before submitting to maximize the chance that your experiment runs
% successfully.
r.debug();


% Marker Size Helper Function to scale size of markers for robots with figure window
% Input: robotarium class instance
function marker_size = determine_robot_marker_size(robotarium_instance)

% Get the size of the robotarium figure window in pixels
curunits = get(robotarium_instance.figure_handle, 'Units');
set(robotarium_instance.figure_handle, 'Units', 'Pixels');
cursize = get(robotarium_instance.figure_handle, 'Position');
set(robotarium_instance.figure_handle, 'Units', curunits);

% Determine the ratio of the robot size to the x-axis (the axis are
% normalized so you could do this with y and figure height as well).
robot_ratio = (robotarium_instance.robot_diameter + 0.03)/...
    (robotarium_instance.boundaries(2) - robotarium_instance.boundaries(1));

% Determine the marker size in points so it fits the window. cursize(3) is
% the width of the figure window in pixels. (the axis are
% normalized so you could do this with y and figure height as well).
marker_size = cursize(3) * robot_ratio;

end

% Font Size Helper Function to scale size with figure window
% Input: robotarium instance, desired height of the font in meters
function font_size = determine_font_size(robotarium_instance, font_height_meters)

% Get the size of the robotarium figure window in point units
curunits = get(robotarium_instance.figure_handle, 'Units');
set(robotarium_instance.figure_handle, 'Units', 'Pixels');
cursize = get(robotarium_instance.figure_handle, 'Position');
set(robotarium_instance.figure_handle, 'Units', curunits);

% Determine the ratio of the font height to the y-axis
font_ratio = (font_height_meters)/(robotarium_instance.boundaries(4) -...
    robotarium_instance.boundaries(3));

% Determine the font size in points so it fits the window. cursize(4) is
% the hight of the figure window in points.
font_size = cursize(4) * font_ratio;

end