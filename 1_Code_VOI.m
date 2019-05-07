% ##############################################################################################%
%
% MODEL CODE FOR:
%   "General rules for environmental management to prioritise social-ecological
%   systems research based on a value of information approach"
%
% AUTHORS:
%   Katrina J Davis, Iadine Chadès, Jonathan R Rhodes & Michael Bode
%
% DESCRIPTION:
%   Code calculates VOI for individual parameters (Analysis 1)
%   and for paired parameters (Analysis 2) in all systems
%
% Note: this code should be run before '2_Code_EVPXI'
% Note: all discretisations have been reduced (as below) for tractability:
%           Replication = 10 (here 'o' and 'u')
%           n = 5 (discretisations for parameter of interest, here 'DiscSteps')
%           b = 10 (discretisations other parameters, here 'MaxP')
% ##############################################################################################%


%% ##############################################################################################%
% ANALYSIS 1 - SYSTEM 1

clear all, tic
rng(1,'twister'); % set random number seed and generator

% PARAMETER COMBINATION;
% I = Influence: 0 < I < 1
% q = Willingness to get on board: 0 < q < 1
% r = Growth rate: 0.2 < r < 2
% C = Connectivity: 0 < C < 1
% H = Harvest rate: 0 < H < 0.4
D = [1 1]*0.25; % Intervention impact

% Simulated Outcome: all outcomes of intervention, influence and engagement
P = [1 1; 0 1; 1 0; 0 0];

% Counter for multiple replications
o = 10;
for w = 1:o % Loop to generate multi-dimensional array with result replication
    
    % % These two parameters define the computational intensity of the analysis
    DiscSteps = 5; % Range for each 'parameter of interest' dimension is split into this many values
    MaxP = 10; % We run this many random replicates of the 'other parameters'
    
    % % Pre-set the random parameters that we're searching across
    % Parameter of interest
    LoopValues = linspace(1e-2,1-1e-2,DiscSteps);
    [x,y] = meshgrid(LoopValues,LoopValues); x = x(:); y = y(:);
    
    % Other parameters
    V_Log = repmat(rand(MaxP,10).*repmat([1 1 1 1 1.8 1.8 1 1 0.4 0.4],MaxP,1) ...
        + repmat([0 0 0 0 0.2 0.2 0 0 0 0],MaxP,1),length(x),1);
    
    for ParameterInterest = 1:5
        disp(['Parameter of interest #' num2str(ParameterInterest)])
        for i = 1:DiscSteps.^2
            if mod(i,10) == 0; disp(['Finished ' num2str(i./DiscSteps.^2) '%']); end
            
            for ParameterCombination = 1:MaxP
                
                % Assigns values of 'other parameters' from V_Log matrix
                % We run comparable random parameters for all steps
                V = V_Log(ParameterCombination + MaxP*(i-1),:);
                I = [1 V(1); V(2) 1]; q = [V(3) V(4)]; r = [V(5) V(6)];
                C = [1-V(7) V(7); V(8) 1-V(8)]; H = [V(9) V(10)];
                
                % Choose the parameter that's being targeted
                % Cycles through range of parameter of interest
                if ParameterInterest == 1
                    I  = [nan x(i); y(i) nan];
                elseif ParameterInterest == 2
                    q = [x(i) y(i)];
                elseif ParameterInterest == 3
                    r = [x(i) y(i)]*1.8 + 0.2;
                elseif ParameterInterest == 4
                    C = [1-x(i) x(i); y(i) 1-y(i)];
                elseif ParameterInterest == 5
                    H = [x(i) y(i)].*0.4;
                end
                
                for A = 1:2; % INTERVENTION
                    % Outcomes = [1 1; 0 1; 1 0; 0 0]
                    if A == 1 % Calculate probability of different outcomes if intervene in S1
                        p_P = [q(1)*I(1,2)*q(2); ...
                            0; ...
                            q(1)*(1 - q(2)*I(1,2)); ...
                            1-q(1)];
                    elseif A == 2 % Calculate probability of different outcomes if intervene in S2
                        p_P = [q(2)*I(2,1)*q(1); ...
                            q(2)*(1 - q(1)*I(2,1)); ...
                            0; ...
                            1-q(2)];
                    end
                    % Ecological Model
                    for p = 1:4
                        Change = 1; dt = 0.2; N = [0.5 0.5];  % Initialise populations
                        % After this the most recent equilibrium is used as the starting point
                        while Change > 1e-6 % Run the population dynamics to equilibrium
                            N_old = N;
                            N(1) = min(1e1,max(0,N(1) + dt*(r(1)*N(1)*(1 - N(1))... % Growth rate of the pop
                                - C(1,2)*N(1) + C(2,1)*N(2)... % Connectivity between the populations
                                - H(1)*(1 - D(1)*P(p,1))*N(1)))); % Harvest of N1 by group 1
                            N(2) = min(1e1,max(0,N(2) + dt*(r(2)*N(2)*(1 - N(2))...
                                + C(1,2)*N(1) - C(2,1)*N(2)...
                                - H(2)*(1 - D(2)*P(p,2))*N(2))));
                            Change = sum(abs(N_old - N));
                        end
                        N_star(p,:) = N;
                    end
                    % Performance if we know perfect information (for each action)
                    % We standardise performance against the system for each parameter set separately.
                    SN = sum(N_star,2) - min(sum(N_star,2));
                    SN = SN./max(SN);
                    PI(ParameterCombination + MaxP*(i-1),A) = ...
                        p_P'*SN;
                end
            end
            PI_X(i,:) = mean(PI(MaxP*(i-1)+1:MaxP*i,:));
        end
        Action_under_complete_certainty(ParameterInterest) = mean(max(PI,[],2))
        Action_under_complete_uncertainty(ParameterInterest) = max(mean(PI))
        Action_under_partially_resolved(ParameterInterest) = mean(max(PI_X,[],2))
    end
    Action_Outcomes(:,:,w) = [Action_under_complete_uncertainty; ...
        Action_under_partially_resolved; ...
        Action_under_complete_certainty];
end

save S1_Single_EVPXI Action_Outcomes

%% ##############################################################################################%
% % ANALYSIS 1 - SYSTEM 2

clear all, tic
rng(1,'twister'); % set the random number seed and generator

% PARAMETER COMBINATION;
% I = Influence: 0 < I < 1
% q = Willingness to get on board: 0 < q < 1
% r = Growth rate: 0.2 < r < 2
% C = Connectivity: 0 < C < 1
% H = Harvest rate: 0 < H < 0.4
M = [1 1]*0.25; % Manager's addition rate
D = [1 1]*0.25; % Intervention impact

% Simulated Outcome
P = [1 1 0 1; 1 1 0 0; 1 0 0 0; 0 1 1 1; 0 0 1 1; 0 0 1 0];

% Counter for multiple replications
o = 10;

for w = 1:o % Loop to generate multi-dimensional array with results
    
    % % These two parameters define the computational intensity of the analysis
    DiscSteps = 5; % Range for each 'parameter of interest' dimension is split into this many values
    MaxP = 10; % We run this many random replicates of the 'other parameters'
    
    % % Pre-set the random parameters that we're searching across
    % Parameter of interest
    LoopValues = linspace(1e-2,1-1e-2,DiscSteps);
    [x,y] = meshgrid(LoopValues,LoopValues); x = x(:); y = y(:);
    
    % Other parameters
    V_Log = repmat(rand(MaxP,10).*repmat([1 1 1 1 1.8 1.8 1 1 0.4 0.4],MaxP,1) ...
        + repmat([0 0 0 0 0.2 0.2 0 0 0 0],MaxP,1),length(x),1);
    
    for ParameterInterest = 1:5
        disp(['Parameter of interest #' num2str(ParameterInterest)])
        for i = 1:DiscSteps.^2
            if mod(i,10) == 0; disp(['Finished ' num2str(i./DiscSteps.^2) '%']); end
            
            for ParameterCombination = 1:MaxP
                V = V_Log(ParameterCombination + MaxP*(i-1),:);
                
                I = [1 V(1); V(2) 1]; q = [V(3) V(4)]; r = [V(5) V(6)];
                C = [1-V(7) V(7); V(8) 1-V(8)]; H = [V(9) V(10)];
                
                % Choose the parameter that's being targeted
                % Cycles through values of parameter of interest
                if ParameterInterest == 1
                    I  = [nan x(i); y(i) nan];
                elseif ParameterInterest == 2
                    q = [x(i) y(i)];
                elseif ParameterInterest == 3
                    r = [x(i) y(i)]*1.8 + 0.2;
                elseif ParameterInterest == 4
                    C = [1-x(i) x(i); y(i) 1-y(i)];
                elseif ParameterInterest == 5
                    H = [x(i) y(i)].*0.4;
                end
                
                for A = 1:2; % INTERVENTION
                    if A == 1
                        p_P = [q(1)*I(1,2)*q(2); ...
                            q(1)*( I(1,2)*(1 - q(2)) + 1 - I(1,2)); ...
                            (1 - q(1)); ...
                            0; ...
                            0; ...
                            0];
                    elseif A == 2
                        p_P = [0; ...
                            0; ...
                            0; ...
                            q(2)*I(2,1)*q(1); ...
                            q(2)*(I(2,1)*(1 - q(1)) + 1 - I(2,1)); ...
                            (1 - q(2))];
                    end
                    % Ecological Model
                    for p = 1:6
                        Change = 1; dt = 0.2; N = [0.5 0.5];
                        while Change > 1e-6
                            N_old = N;
                            N(1) = min(1e1,max(0,N(1) + dt*(r(1)*N(1)*(1 - N(1))... % Growth rate of the pop
                                - C(1,2)*N(1) + C(2,1)*N(2)... % Connectivity between the populations
                                + M(1)*P(p,1)*N(1) ... % Supplement of N1 by manager
                                - H(1)*(1 + D(1)*P(p,2))*N(1)))); % Harvest of N1 by S1
                            N(2) = min(1e1,max(0,N(2) + dt*(r(2)*N(2)*(1 - N(2))...
                                + C(1,2)*N(1) - C(2,1)*N(2)...
                                + M(2)*P(p,3)*N(2)...
                                - H(2)*(1 + D(2)*P(p,4))*N(2))));
                            Change = sum(abs(N_old - N));
                        end
                        N_star(p,:) = N;
                    end
                    % Performance if we know the perfect information (for each action)
                    % We standardise performance against the system for each parameter set separately
                    SN = sum(N_star,2) - min(sum(N_star,2));
                    SN = SN./max(SN);
                    PI(ParameterCombination + MaxP*(i-1),A) = ...
                        p_P'*SN;
                end
            end
            PI_X(i,:) = mean(PI(MaxP*(i-1)+1:MaxP*i,:));
        end
        Action_under_complete_certainty(ParameterInterest) = mean(max(PI,[],2))
        Action_under_complete_uncertainty(ParameterInterest) = max(mean(PI))
        Action_under_partially_resolved(ParameterInterest) = mean(max(PI_X,[],2))
    end
    Action_Outcomes(:,:,w) = [Action_under_complete_uncertainty; ...
        Action_under_partially_resolved; ...
        Action_under_complete_certainty];
end

save S2_Single_EVPXI Action_Outcomes

%% ##############################################################################################%
% % ANALYSIS 1 - SYSTEM 3

clear all, tic
rng(1,'twister'); % set the random number seed and generator

% PARAMETER COMBINATION;
% I = Influence: 0 < I < 1
% q = Willingness to get on board: 0 < q < 1
% r = Growth rate: 0.2 < r < 2
% C = Connectivity: 0 < C < 1
% H = Harvest rate: 0 < H < 0.4
M = [1 1]*0.25; % Manager's addition rate
D = [1 1]*0.25; % Intervention impact

% Simulated Outcome
P = [1 1 0 1; 1 1 0 0; 1 0 0 0; 0 1 1 1; 0 0 1 1; 0 0 1 0];

% Counter for multiple replications
o = 10;

for w = 1:o % Loop to generate multi-dimensional array with results
    % These two parameters define the computational intensity of the analysis
    DiscSteps = 5 % Range for each 'parameter of interest' dimension is split into this many values
    MaxP = 10 % We run this many random replicates of the 'other parameters'
    
    % % Pre-set the random parameters that we're searching across
    % Parameter of interest
    LoopValues = linspace(1e-2,1-1e-2,DiscSteps);
    [x,y] = meshgrid(LoopValues,LoopValues); x = x(:); y = y(:);
    
    % Other parameters
    V_Log = repmat(rand(MaxP,10).*repmat([1 1 1 1 1.8 1.8 1 1 0.4 0.4],MaxP,1) ...
        + repmat([0 0 0 0 0.2 0.2 0 0 0 0],MaxP,1),length(x),1);
    
    for ParameterInterest = 1:5
        disp(['Parameter of interest #' num2str(ParameterInterest)])
        for i = 1:DiscSteps.^2
            if mod(i,10) == 0; disp(['Finished ' num2str(i./DiscSteps.^2) '%']); end
            
            for ParameterCombination = 1:MaxP
                
                % Assigns values of 'other parameters' from V_Log matrix
                % We run comparable random parameters for all steps
                V = V_Log(ParameterCombination + MaxP*(i-1),:);
                I = [1 V(1); V(2) 1]; q = [V(3) V(4)]; r = [V(5) V(6)];
                C = [1-V(7) V(7); 1-V(8) V(8)]; H = [V(9) V(10)];
                %                 M = [V(11) V(12)];
                
                % Choose the parameter that's being targeted
                % Cycles through values of parameter of interest
                if ParameterInterest == 1
                    I  = [nan x(i); y(i) nan];
                elseif ParameterInterest == 2
                    q = [x(i) y(i)];
                elseif ParameterInterest == 3
                    r = [x(i) y(i)]*1.8 + 0.2;
                elseif ParameterInterest == 4
                    C = [1-x(i) x(i); 1-y(i) y(i)];
                elseif ParameterInterest == 5
                    H = [x(i) y(i)].*0.4;
                end
                
                for A = 1:2; % INTERVENTION
                    if A == 1
                        p_P = [q(1)*I(1,2)*q(2); ...
                            q(1)*( I(1,2)*(1 - q(2)) + 1 - I(1,2)); ...
                            (1 - q(1)); ...
                            0; ...
                            0; ...
                            0];
                    elseif A == 2
                        p_P = [0; ...
                            0; ...
                            0; ...
                            q(2)*I(2,1)*q(1); ...
                            q(2)*(I(2,1)*(1 - q(1)) + 1 - I(2,1)); ...
                            (1 - q(2))];
                    end
                    
                    % Ecological Model
                    for p = 1:6
                        Change = 1; dt = 0.2; N = [0.5 0.5];
                        while Change > 1e-6
                            N_old = N;
                            N(1) = min(1e1,max(0,N(1) + dt*(r(1)*N(1)*(1 - N(1))... % Growth rate of the pop
                                - C(1,2)*N(1) + C(2,1)*N(2)... % Connectivity between the populations
                                - M(1)*P(p,1)*N(1)... % Harvest of N1 by manager
                                - H(1)*(1 - D(1)*P(p,2))*N(1)))); % Harvest of N1 by farmer 1
                            N(2) = min(1e1,max(0,N(2) + dt*(r(2)*N(2)*(1 - N(2))...
                                + C(1,2)*N(1) - C(2,1)*N(2)...
                                - M(2)*P(p,3)*N(2)...
                                - H(2)*(1 - D(2)*P(p,4))*N(2))));
                            Change = sum(abs(N_old - N));
                        end
                        
                        % GINI Coefficient
                        a = 0.5; % GINI - cumulative proportion of pop from S1
                        b = min(N)/(N(1)+ N(2)); % GINI - cumulative share
                        L = 0.5 *(-a+b+1); % GINI - area under the Lorenz curve
                        E = 0.5 - L; % GINI - Area under line of equality and above Lorenz curve
                        Gini = E/(E+L); % GINI - gini coefficient
                        %                         N_star(p,:) = N; % Old N_star
                        N_starG(p,:) = (N(1) + N(2))*(1-Gini);
                        
                    end
                    % Performance if we know the perfect information (for each action)
                    % We standardise performance against the system for each parameter set separately.
                    SN = max(N_starG)-N_starG;
                    SN = SN./max(SN);
                    PI(ParameterCombination + MaxP*(i-1),A) = p_P'*SN;
                end
            end
            PI_X(i,:) = mean(PI(MaxP*(i-1)+1:MaxP*i,:));
        end
        Action_under_complete_certainty(ParameterInterest) = mean(max(PI,[],2))
        Action_under_complete_uncertainty(ParameterInterest) = max(mean(PI))
        Action_under_partially_resolved(ParameterInterest) = mean(max(PI_X,[],2))
    end
    Action_Outcomes(:,:,w) = [Action_under_complete_uncertainty; ...
        Action_under_partially_resolved; ...
        Action_under_complete_certainty];
end
save S3_Single_EVPXI Action_Outcomes


%% ##############################################################################################%
% % ANALYSIS 1 - SYSTEM 4

clear all, tic
rng(1,'twister'); % set the random number seed and generator

% PARAMETER COMBINATION;
% I = Influence: 0 < I < 1
% q = Willingness to get on board: 0 < q < 1
% r = Growth rate: 0.2 < r < 2
% C = Connectivity: 0 < C < 1
% H = Harvest rate: 0 < H < 0.4
D = [1 1]*0.25; % Intervention impact

% Simulated Outcome
P = [1 1; 0 1; 1 0; 0 0];

% Counter for multiple replications
o = 10;

for w = 1:o % Loop to generate multi-dimensional array with results
    
    % These two parameters define the computational intensity of the analysis
    DiscSteps = 5; % Range for each 'parameter of interest' dimension is split into this many values
    MaxP = 10; % We run this many random replicates of the 'other parameters'
    
    % % Pre-set the random parameters that we're searching across
    % Parameter of interest
    LoopValues = linspace(1e-2,1-1e-2,DiscSteps);
    [x,y] = meshgrid(LoopValues,LoopValues); x = x(:); y = y(:);
    
    % Other parameters
    V_Log = repmat(rand(MaxP,10).*repmat([1 1 1 1 1.8 1.8 1 1 0.4 0.4],MaxP,1) ...
        + repmat([0 0 0 0 0.2 0.2 0 0 0 0],MaxP,1),length(x),1);
    
    for ParameterInterest = 1:5
        disp(['Parameter of interest #' num2str(ParameterInterest)])
        for i = 1:DiscSteps.^2
            if mod(i,10) == 0; disp(['Finished ' num2str(i./DiscSteps.^2) '%']); end
            
            for ParameterCombination = 1:MaxP
                
                % Assigns values of 'other parameters' from V_Log matrix
                % We run comparable random parameters for all steps
                V = V_Log(ParameterCombination + MaxP*(i-1),:);
                I = [1 V(1); V(2) 1]; q = [V(3) V(4)]; r = [V(5) V(6)];
                C = [1-V(7) V(7); 1-V(8) V(8)]; H = [V(9) V(10)];
                
                % Choose the parameter that's being targeted here
                % Cycles through values of parameter of interest
                if ParameterInterest == 1
                    I  = [nan x(i); y(i) nan];
                elseif ParameterInterest == 2
                    q = [x(i) y(i)];
                elseif ParameterInterest == 3
                    r = [x(i) y(i)]*1.8 + 0.2;
                elseif ParameterInterest == 4
                    C = [1-x(i) x(i); 1-y(i) y(i)];
                elseif ParameterInterest == 5
                    H = [x(i) y(i)].*0.4;
                end
                
                for A = 1:2; % INTERVENTION
                    if A == 1
                        p_P = [q(1)*q(2)*I(1,2); ...
                            0; ...
                            q(1)*(1 - q(2)*I(1,2)); ...
                            1-q(1)];
                    elseif A == 2
                        p_P = [q(2)*q(1)*I(2,1); ...
                            q(2)*(1 - q(1)*I(2,1)); ...
                            0; ...
                            1-q(2)];
                    end
                    % Ecological Model
                    for p = 1:4
                        Change = 1; N = [0.5 0.5]; dt = 0.2;
                        while Change > 1e-6
                            N_old = N;
                            N(1) = min(1e1,max(0,N(1) + dt*(r(1)*N(1)*(1 - N(1))... % Growth rate of the pop
                                - C(1,2)*N(1) + C(2,1)*N(2)... % Connectivity between the populations
                                - H(1)*(1 - D(1)*P(p,1))*N(1))));  % Harvest of N1 by S1
                            N(2) = min(1e1,max(0,N(2) + dt*(r(2)*N(2)*(1 - N(2))...
                                + C(1,2)*N(1) - C(2,1)*N(2)...
                                - H(2)*(1 - D(2)*P(p,2))*N(2))));
                            Change = sum(abs(N_old - N));
                        end
                        % GINI Coefficient
                        a = 0.5; % cumulative proportion of pop from S1
                        b = min(N)/(N(1)+ N(2)); % cumulative share
                        L = 0.5*(-a + b + 1); % area under the Lorenz curve
                        E = 0.5 - L; % area under line of equality and above Lorenz curve
                        Gini = E/(E+L); % gini coefficient
                        %                         N_star(p,:) = N; % Old N_star
                        N_starG(p,:) = (N(1) + N(2))*(1-Gini);
                    end
                    
                    % Performance if we know the perfect information (for each action)
                    % We standardise performance against the system, defined by the
                    % parameters.
                    SN = sum(N_starG,2) - min(sum(N_starG,2));
                    SN = SN./max(SN);
                    PI(ParameterCombination + MaxP*(i-1),A) = ...
                        p_P'*SN;
                end
            end
            PI_X(i,:) = mean(PI(MaxP*(i-1)+1:MaxP*i,:));
        end
        Action_under_complete_certainty(ParameterInterest) = mean(max(PI,[],2));
        Action_under_complete_uncertainty(ParameterInterest) = max(mean(PI));
        Action_under_partially_resolved(ParameterInterest) = mean(max(PI_X,[],2));
    end
    
    Action_Outcomes(:,:,w) = [Action_under_complete_uncertainty; ...
        Action_under_partially_resolved; ...
        Action_under_complete_certainty];
end
save S4_Single_EVPXI Action_Outcomes


%% ##############################################################################################%
% % ANALYSIS 2 - SYSTEM 1

clear all, tic
rng(1,'twister'); % set the random number seed and generator

% PARAMETER COMBINATION;
% I = Influence: 0 < I < 1
% q = Willingness to get on board: 0 < q < 1
% r = Growth rate: 0.2 < r < 2 
% C = Connectivity: 0 < C < 1
% H = Harvest rate: 0 < H < 0.4 
D = [1 1]*0.25; % Intervention impact

% Simulated Outcome
P = [1 1; 0 1; 1 0; 0 0];

% Counter for multiple replications
u = 10;

for t = 1:u % Loop to generate multi-dimensional array with with result replication
   
    % % These two parameters define the computational intenstity of the analysis
    DiscSteps = 5; % Range for each 'parameter of interest' dimension is split into this many values
    MaxP = 10; % We run this many random replicates of the 'other parameters'
    
    % % Pre-set the random parameters that we're searching across.
    %   There are four, since each parameter is really two parameters (e.g., q1 & q2)
    % Parameter of interest
    LoopValues = linspace(1e-2,1-1e-2,DiscSteps);
    [w,x,y,z] = ndgrid(LoopValues,LoopValues,LoopValues,LoopValues);
    x = x(:); y = y(:); w = w(:); z = z(:);
    
    % Other parameters
    UseSame = 1;
    if UseSame == 1
        V_Log = repmat(rand(MaxP,10).*repmat([1 1 1 1 1.8 1.8 1 1 0.4 0.4],MaxP,1) ...
            + repmat([0 0 0 0 0.2 0.2 0 0 0 0],MaxP,1),length(x),1); 
    end
    
    for ParameterInterest = 1:3
        disp(['Parameter of interest #' num2str(ParameterInterest)])
        for i = 1:DiscSteps.^4 
            if mod(i,100) == 0; disp(['Finished ' num2str(i./DiscSteps.^4) '%']); end
            
            % Parameters combinations
            for ParameterCombination = 1:MaxP
                if UseSame == 1
                    V = V_Log(ParameterCombination + MaxP*(i-1),:);
                else
                    V = rand(1,10).*[1 1 1 1 1.8 1.8 1 1 0.4 0.4] + [0 0 0 0 0.2 0.2 0 0 0 0];
                end
                I = [1 V(1); V(2) 1]; q = [V(3) V(4)]; r = [V(5) V(6)];
                C = [1-V(7) V(7); V(8) 1-V(8)]; H = [V(9) V(10)];
                
                % Choose parameters that are being targeted
                if ParameterInterest == 1 % Social parameters
                    I = [nan w(i); x(i) nan];
                    q = [y(i) z(i)];
                elseif ParameterInterest == 2 % Ecological parameters
                    r = [w(i) x(i)]*1.8 + 0.2;
                    C = [1-y(i) y(i); z(i) 1-z(i)];
                elseif ParameterInterest == 3 % Socio-ecological linkages
                    H = [w(i) y(i)].*0.4;
                end
                
                for A = 1:2; % INTERVENTION
                    if A == 1 
                        p_P = [q(1)*q(2)*I(1,2); ...
                            0; ...
                            q(1)*(1 - q(2)*I(1,2)); ...
                            1-q(1)];
                    elseif A == 2
                        p_P = [q(2)*q(1)*I(2,1); ...
                            q(2)*(1 - q(1)*I(2,1)); ...
                            0; ...
                            1-q(2)];
                    end
                    % Ecological Model
                    for p = 1:4
                        Change = 1; dt = 0.2; N = [0.5 0.5];
                        while Change > 1e-6
                            N_old = N;
                            N(1) = min(1e1,max(0,N(1) + dt*(r(1)*N(1)*(1 - N(1))... % growth rate of the pop
                                - C(1,2)*N(1) + C(2,1)*N(2)... % connectivity between the populations
                                - H(1)*(1 - D(1)*P(p,1))*N(1)))); % harvest of N1 by group 1
                            N(2) = min(1e1,max(0,N(2) + dt*(r(2)*N(2)*(1 - N(2))...
                                + C(1,2)*N(1) - C(2,1)*N(2)...
                                - H(2)*(1 - D(2)*P(p,2))*N(2))));
                            Change = sum(abs(N_old - N));
                        end
                        N_star(p,:) = N;
                    end
                    % Performance if we know the perfect information (for each action)
                    % We standardise performance against the system for each parameter set separately.
                    SN = sum(N_star,2) - min(sum(N_star,2));
                    SN = SN./max(SN);
                    PI(ParameterCombination + MaxP*(i-1),A) = ...
                        p_P'*SN;
                end
            end
            PI_X(i,:) = mean(PI(MaxP*(i-1)+1:MaxP*i,:));
        end
        Action_under_complete_certainty(ParameterInterest) = mean(max(PI,[],2))
        Action_under_complete_uncertainty(ParameterInterest) = max(mean(PI))
        Action_under_partially_resolved(ParameterInterest) = mean(max(PI_X,[],2))
    end
    
    Action_DblOutcomes(:,:,t) = [Action_under_complete_uncertainty; ...
        Action_under_partially_resolved; ...
        Action_under_complete_certainty]
end
save S1_Pairwise_EVPXI Action_DblOutcomes

%% ##############################################################################################%
% % ANALYSIS 2 - SYSTEM 2
% ##############################################################################################%

clear all, tic
rng(1,'twister'); % set the random number seed and generator

% PARAMETER COMBINATION;
% I = Influence: 0 < I < 1
% q = Willingness to get on board: 0 < q < 1
% r = Growth rate: 0.2 < r < 2 
% C = Connectivity: 0 < C < 1
% H = Harvest rate: 0 < H < 0.4 
M = [1 1]*0.25; % Manager's addition rate
D = [1 1]*0.25; % Intervention impact

% Simulated Outcome
P = [1 1 0 1; 1 1 0 0; 1 0 0 0; 0 1 1 1; 0 0 1 1; 0 0 1 0];

% Counter for multiple replications
u = 10

for t = 1:u % Loop to generate multi-dimensional array with with result replication
% These two parameters define the computational intenstity of the analysis
    DiscSteps = 5 %20; % Range for each 'parameter of interest' dimension is split into this many values [20]
    MaxP = 10 %250; % We run this many random replicates of the 'other parameters' [250]
    
    % % Pre-set the random parameters that we're searching across.
    %   There are four, since each parameter is really two parameters (e.g., q1 & q2)
    % Parameter of interest
    LoopValues = linspace(1e-2,1-1e-2,DiscSteps);
    [w,x,y,z] = ndgrid(LoopValues,LoopValues,LoopValues,LoopValues);
    x = x(:); y = y(:); w = w(:); z = z(:);
    
    % Other parameters
    UseSame = 1;
    if UseSame == 1
        V_Log = repmat(rand(MaxP,10).*repmat([1 1 1 1 1.8 1.8 1 1 0.4 0.4],MaxP,1) ...
            + repmat([0 0 0 0 0.2 0.2 0 0 0 0],MaxP,1),length(x),1);
    end
    
    for ParameterInterest = 1:3
        disp(['Parameter of interest #' num2str(ParameterInterest)])
        for i = 1:DiscSteps.^4
            if mod(i,10) == 0; disp(['Finished ' num2str(i./DiscSteps.^4) '%']); end
            
            for ParameterCombination = 1:MaxP
                if UseSame == 1
                    V = V_Log(ParameterCombination + MaxP*(i-1),:);
                else
                    V = rand(1,10).*[1 1 1 1 1.8 1.8 1 1 0.4 0.4] + [0 0 0 0 0.2 0.2 0 0 0 0];
                end
                I = [1 V(1); V(2) 1];           q = [V(3) V(4)];       r = [V(5) V(6)];
                C = [1-V(7) V(7); 1-V(8) V(8)]; H = [V(9) V(10)];
                
                % Choose the parameter that's being targeted
                % Cycles through values of parameter of interest
                if ParameterInterest == 1 % Social parameters
                    I = [nan w(i); x(i) nan];
                    q = [y(i) z(i)];
                elseif ParameterInterest == 2 % Ecological parameters
                    r = [w(i) x(i)]*1.8 + 0.2; 
                    C = [1-y(i) y(i); 1-z(i) z(i)];
                elseif ParameterInterest == 3 % Socio-ecological linkages
                    H = [w(i) x(i)].*0.4;
                end
                
                for A = 1:2; % Intervention
                    if A == 1 
                        p_P = [q(1)*I(1,2)*q(2); ...
                            q(1)*( I(1,2)*(1 - q(2)) + 1 - I(1,2)); ...
                            (1 - q(1)); ...
                            0; ...
                            0; ...
                            0];
                    elseif A == 2
                        p_P = [0; ...
                            0; ...
                            0; ...
                            q(2)*I(2,1)*q(1); ...
                            q(2)*(I(2,1)*(1 - q(1)) + 1 - I(2,1)); ...
                            (1 - q(2))];
                    end
                    %% Ecological Model
                    for p = 1:6
                        Change = 1; dt = 0.2; N = [0.5 0.5];
                        while Change > 1e-6
                            N_old = N;
                            N(1) = min(1e1,max(0,N(1) + dt*(r(1)*N(1)*(1 - N(1))... % Growth rate of the pop
                                - C(1,2)*N(1) + C(2,1)*N(2)... % Connectivity between the populations
                                + M(1)*P(p,1)*N(1) ... % Supplement of N1 by manager
                                - H(1)*(1 + D(1)*P(p,2))*N(1)))); % Harvest of N1 by S1
                            N(2) = min(1e1,max(0,N(2) + dt*(r(2)*N(2)*(1 - N(2))...
                                + C(1,2)*N(1) - C(2,1)*N(2)...
                                + M(2)*P(p,3)*N(2)...
                                - H(2)*(1 + D(2)*P(p,4))*N(2))));
                            Change = sum(abs(N_old - N));
                        end
                        N_star(p,:) = N;
                    end
                    % Performance if we know the perfect information (for each action)
                    % We standardise performance against the system for each parameter set separately.
                    SN = sum(N_star,2) - min(sum(N_star,2));
                    SN = SN./max(SN);
                    PI(ParameterCombination + MaxP*(i-1),A) = ...
                        p_P'*SN;
                end
            end
            PI_X(i,:) = mean(PI(MaxP*(i-1)+1:MaxP*i,:)); 
        end
        Action_under_complete_certainty(ParameterInterest) = mean(max(PI,[],2))
        Action_under_complete_uncertainty(ParameterInterest) = max(mean(PI))
        Action_under_partially_resolved(ParameterInterest) = mean(max(PI_X,[],2))
    end
    
    Action_DblOutcomes(:,:,t) = [Action_under_complete_uncertainty; ...
        Action_under_partially_resolved; ...
        Action_under_complete_certainty];    
end
save S2_Pairwise_EVPXI Action_DblOutcomes


%% ##############################################################################################%
% % ANALYSIS 2 - SYSTEM 3

clear all, tic
rng(1,'twister'); % set the random number seed and generator

% PARAMETER COMBINATION;
    % I = Influence: 0 < I < 1
    % q = Willingness to get on board: 0 < q < 1
    % r = Growth rate: 0.2 < r < 2 
    % C = Connectivity: 0 < C < 1
    % H = Harvest rate: 0 < H < 0.4 
    M = [1 1]*0.25; % Manager's addition rate 
    D = [1 1]*0.25; % Intervention impact

% Simulated outcomes
P = [1 1 0 1; 1 1 0 0; 1 0 0 0; 0 1 1 1; 0 0 1 1; 0 0 1 0];

% Counter for multiple replications
u = 10 %100; % No of loops

for t = 1:u % loop to generate multi-dimensional array with results
    % These two parameters define the computational intenstity of the analysis
    DiscSteps = 5; % Range for each 'parameter of interest' dimension is split into this many values
    MaxP = 10; % We run this many random replicates of the 'other parameters'
    
    % % Pre-set the random parameters that we're searching across.
    %   There are four, since each parameter is really two parameters (e.g., q1 & q2)
    % Parameter of interest
    LoopValues = linspace(1e-2,1-1e-2,DiscSteps);
    [w,x,y,z] = ndgrid(LoopValues,LoopValues,LoopValues,LoopValues);
    x = x(:); y = y(:); w = w(:); z = z(:);
    
    % Other parameters
    UseSame = 1;
    if UseSame == 1
        V_Log = repmat(rand(MaxP,10).*repmat([1 1 1 1 1.8 1.8 1 1 0.4 0.4],MaxP,1) ...
            + repmat([0 0 0 0 0.2 0.2 0 0 0 0],MaxP,1),length(x),1); 
    end
    
    for ParameterInterest = 1:3 
        disp(['Parameter of interest #' num2str(ParameterInterest)])
        for i = 1:DiscSteps.^4
            if mod(i,10) == 0; disp(['Finished ' num2str(i./DiscSteps.^4) '%']); end
            
            
            for ParameterCombination = 1:MaxP 
                if UseSame == 1
                    V = V_Log(ParameterCombination + MaxP*(i-1),:);
                else
                    V = rand(1,10).*[1 1 1 1 1.8 1.8 1 1 0.4 0.4] + [0 0 0 0 0.2 0.2 0 0 0 0];
                end
                I = [1 V(1); V(2) 1];           q = [V(3) V(4)];       r = [V(5) V(6)];
                C = [1-V(7) V(7); 1-V(8) V(8)]; H = [V(9); V(10)];
                
                % Choose the parameter that's being targeted
                % Cycles through values of parameter of interest
                if ParameterInterest == 1 % Social parameters
                    I = [nan w(i); x(i) nan];
                    q = [y(i) z(i)];
                elseif ParameterInterest == 2 % Ecological parameters
                    r = [w(i) x(i)]*1.8 + 0.2;
                    C = [1-y(i) y(i); 1-z(i) z(i)];
                elseif ParameterInterest == 3 % Socio-ecological linkages
                    H = [w(i) x(i)].*0.4;
                end
                
                for A = 1:2; % Intervention
                    if A == 1 
                        p_P = [q(1)*I(1,2)*q(2); ...
                            q(1)*( I(1,2)*(1 - q(2)) + 1 - I(1,2)); ...
                            (1 - q(1)); ...
                            0; ...
                            0; ...
                            0];
                    elseif A == 2
                        p_P = [0; ...
                            0; ...
                            0; ...
                            q(2)*I(2,1)*q(1); ...
                            q(2)*(I(2,1)*(1 - q(1)) + 1 - I(2,1)); ...
                            (1 - q(2))];
                    end
                    % Ecological Model
                    for p = 1:6
                        Change = 1; dt = 0.2; N = [0.5 0.5];
                        while Change > 1e-6
                            N_old = N;
                            N(1) = min(1e1,max(0,N(1) + dt*(r(1)*N(1)*(1 - N(1))... % Growth rate of the pop
                                - C(1,2)*N(1) + C(2,1)*N(2)... % Connectivity between the populations
                                - M(1)*P(p,1)*N(1) ... % Harvest by manager
                                - H(1)*(1 - D(1)*P(p,2))*N(1)))); % Harvest of N1 by S1
                            N(2) = min(1e1,max(0,N(2) + dt*(r(2)*N(2)*(1 - N(2))...
                                + C(1,2)*N(1) - C(2,1)*N(2)...
                                - M(2)*P(p,3)*N(2)... 
                                - H(2)*(1 - D(2)*P(p,4))*N(2))));
                            Change = sum(abs(N_old - N));
                        end
                        % GINI Coefficient 
                        a = 0.5; % cumulative proportion of pop from S1
                        b = min(N)/(N(1)+ N(2)); % cumulative share
                        L = 0.5 *(-a+b+1); % area under the Lorenz curve
                        E = 0.5 - L; % area under line of equality and above Lorenz curve
                        Gini = E/(E+L); % gini coefficient
                        N_starG(p,:) = (N(1) + N(2))*(1-Gini);
                    end
                    % Performance if we know the perfect information (for each action)
                    % We standardise performance against the system for each parameter set separately.
                    SN = max(N_starG)- N_starG; 
                    SN = SN./max(SN); 
                    PI(ParameterCombination + MaxP*(i-1),A) = p_P'*SN; 
                    
                end
            end
            PI_X(i,:) = mean(PI(MaxP*(i-1)+1:MaxP*i,:));
        end
        Action_under_complete_certainty(ParameterInterest) = mean(max(PI,[],2));
        Action_under_complete_uncertainty(ParameterInterest) = max(mean(PI));
        Action_under_partially_resolved(ParameterInterest) = mean(max(PI_X,[],2));
    end
    
    Action_DblOutcomes(:,:,t) = [Action_under_complete_uncertainty; ...
        Action_under_partially_resolved; ...
        Action_under_complete_certainty];
end
save S3_Pairwise_EVPXI Action_DblOutcomes

%% ##############################################################################################%
% % ANALYSIS 2 - SYSTEM 4

clear all, tic
rng(1,'twister'); % set the random number seed and generator

% PARAMETER COMBINATION;
% I = Influence: 0 < I < 1
% q = Willingness to get on board: 0 < q < 1
% r = Growth rate: 0.2 < r < 2
% C = Connectivity: 0 < C < 1
% H = Harvest rate: 0 < H < 0.4
D = [1 1]*0.25; % Intervention impact

% Simlulated outcomes 
P = [1 1; 0 1; 1 0; 0 0];

% Counter for multiple replication 
u = 10 

for t = 1:u % Loop to generate multi-dimensional array with results
    
    % These two parameters define the computational intenstity of the analysis
    DiscSteps = 5; % Range for each 'parameter of interest' dimension is split into this many values
    MaxP = 10; % We run this many random replicates of the 'other parameters'
    
    % % Pre-set the random parameters that we're searching across.
    %   There are four, since each parameter is really two parameters (e.g., q1 & q2)
    % Parameter of interest
    LoopValues = linspace(1e-2,1-1e-2,DiscSteps);
    [w,x,y,z] = ndgrid(LoopValues,LoopValues,LoopValues,LoopValues);
    x = x(:); y = y(:); w = w(:); z = z(:);
    
    % Other parameters
    UseSame = 1;
    if UseSame == 1
        V_Log = repmat(rand(MaxP,10).*repmat([1 1 1 1 1.8 1.8 1 1 0.4 0.4],MaxP,1) ...
            + repmat([0 0 0 0 0.2 0.2 0 0 0 0],MaxP,1),length(x),1);
    end
    
    for ParameterInterest = 1:3 
        disp(['Parameter of interest #' num2str(ParameterInterest)])
        for i = 1:DiscSteps.^4
            if mod(i,100) == 0; disp(['Finished ' num2str(i./DiscSteps.^4) '%']); end
            
            for ParameterCombination = 1:MaxP
                if UseSame == 1
                    V = V_Log(ParameterCombination + MaxP*(i-1),:);
                else
                    V = rand(1,10).*[1 1 1 1 1.8 1.8 1 1 0.4 0.4] + [0 0 0 0 0.2 0.2 0 0 0 0];
                end
                I = [1 V(1); V(2) 1]; q = [V(3) V(4)]; r = [V(5) V(6)];
                C = [1-V(7) V(7); 1-V(8) V(8)]; H = [V(9) V(10)];
                
                % Choose parameters that are being targeted
                if ParameterInterest == 1 % Social parameters
                    I = [nan w(i); x(i) nan];
                    q = [y(i) z(i)];
                elseif ParameterInterest == 2 % Ecological parameters
                    r = [w(i) x(i)]*1.8 + 0.2;
                    C = [1-y(i) y(i); 1-z(i) z(i)];
                elseif ParameterInterest == 3 % Socio-ecological linkages
                    H = [w(i) y(i)].*0.4;
                end
                
                for A = 1:2; % INTERVENTION
                    if A == 1
                        p_P = [q(1)*q(2)*I(1,2); ...
                            0; ...
                            q(1)*(1 - q(2)*I(1,2)); ...
                            1-q(1)];
                    elseif A == 2
                        p_P = [q(2)*q(1)*I(2,1); ...
                            q(2)*(1 - q(1)*I(2,1)); ...
                            0; ...
                            1-q(2)];
                    end
                    % Ecological Model
                    for p = 1:4
                        Change = 1; dt = 0.2; N = [0.5 0.5];
                        while Change > 1e-6
                            N_old = N;
                            N(1) = min(1e1,max(0,N(1) + dt*(r(1)*N(1)*(1 - N(1))... % Growth rate of the pop
                                - C(1,2)*N(1) + C(2,1)*N(2)... % Connectivity between the populations
                                - H(1)*(1 - D(1)*P(p,1))*N(1)))); % Harvest of N1 by S1
                            N(2) = min(1e1,max(0,N(2) + dt*(r(2)*N(2)*(1 - N(2))...
                                + C(1,2)*N(1) - C(2,1)*N(2)...
                                - H(2)*(1 - D(2)*P(p,2))*N(2))));
                            Change = sum(abs(N_old - N));
                        end
                        % GINI Coefficient
                        a = 0.5; % cumulative proportion of pop from S1
                        b = min(N)/(N(1)+ N(2)); % cumulative share
                        L = 0.5*(-a + b + 1); % area under the Lorenz curve
                        E = 0.5 - L; % area under line of equality and above Lorenz curve
                        Gini = E/(E+L); % gini coefficient
                        N_starG(p,:) = (N(1) + N(2))*(1-Gini);
                     end
                    %% Performance if we know the perfect information (for each action)
                    % We standardise performance against the system for each parameter set separately.
                    SN = sum(N_starG,2) - min(sum(N_starG,2));
                    SN = SN./max(SN);
                    PI(ParameterCombination + MaxP*(i-1),A) = ...
                        p_P'*SN;
                end
            end
            PI_X(i,:) = mean(PI(MaxP*(i-1)+1:MaxP*i,:)); 
        end
        Action_under_complete_certainty(ParameterInterest) = mean(max(PI,[],2));
        Action_under_complete_uncertainty(ParameterInterest) = max(mean(PI));
        Action_under_partially_resolved(ParameterInterest) = mean(max(PI_X,[],2));
    end
    Action_DblOutcomes(:,:,t) = [Action_under_complete_uncertainty; ...
        Action_under_partially_resolved; ...
        Action_under_complete_certainty]
end
save S4_Pairwise_EVPXI Action_DblOutcomes



