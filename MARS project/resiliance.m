%% Consensus with a static, undirected topology & F malicious agents
% MECH 6V29: MARS final project

% Demonstrates a working consensus algorithm where a given number of robots
% will experience failure & remain stationary for a undetermined amount of
% time
clear;clc;

N = 12;
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true);


%% Experiment constants

F = 1;                          % number of malicious agents
G = 1;                          % how many 
iterations = 1000;              % number of iterations the experiment will run over
mal_r = randperm(N,F);          % list of which robots are malicious

% try out different topologies (make sure to change algorithm section
% accordingly)
% L = cycleGL(N);     % constant cycle graph
L = completeGL(N);  % constant complete graph


%% Grab tools we need to convert from single-integrator to unicycle dynamics

% Gain for the diffeomorphism transformation between single-integrator and
% unicycle dynamics
[~, uni_to_si_states] = create_si_to_uni_mapping();
si_to_uni_dyn = create_si_to_uni_dynamics_with_backwards_motion();

uni_barrier_cert_boundary = create_uni_barrier_certificate_with_boundary();

% Initialize velocity vector for agents.  Each agent expects a 2 x 1
% velocity vector containing the linear and angular velocity, respectively.
dxi = zeros(2, N);

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
        
        % do not update velocity if robot is malicious
        if malicious
            continue
        end

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