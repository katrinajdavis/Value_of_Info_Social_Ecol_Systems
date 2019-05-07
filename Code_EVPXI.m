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
%   Code brings together data and calculcates EVPXI for individual parameters (Analysis 1)
%   and for paired parameters (Analysis 2) in all systems
%
% Note: this code should be run after 'Code_VOI'
% ##############################################################################################%

% % PARAMETERS
% 1.   I = Influence
% 2.   q = Willingness to get on board
% 3.   r = Growth rate
% 4.   C = Connectivity
% 5.   H = Harvest rate
% 6.   Social parameters: I, q
% 7.   Ecological parameters: r, C
% 8.   Socio-ecological linkages: H

% % EVPXI RESULTS KEY
% (1) Action_under_complete_uncertainty
% (2) Action_under_partially_resolved
% (3) Action_under_complete_uncertainty

%% System 1

clear all
% Load data
load S1_Single_EVPXI Action_Outcomes
load S1_Pairwise_EVPXI Action_DblOutcomes

% SINGLE PARAMETERS - EVPXI
Performance_single = squeeze((Action_Outcomes(2,:,:) - Action_Outcomes(1,:,:))./...
    (Action_Outcomes(3,:,:) - Action_Outcomes(1,:,:)))';
% Quantiles give us the range we display (max, min and average)
Q_single = quantile(Performance_single,[0.025 0.5 0.975]);

% COMBINATIONS OF PARAMETERS - EVPXI
Performance_double = squeeze((Action_DblOutcomes(2,:,:) - Action_DblOutcomes(1,:,:))./...
    (Action_DblOutcomes(3,:,:) - Action_DblOutcomes(1,:,:)))';
% Quantiles give us the range we display (max, min and average)
Q_double = quantile(Performance_double,[0.025 0.5 0.975]);

% Combine single parameter and combinations of parameters
PA = [Q_single Q_double];

% Write to csv file
csvwrite('S1.csv', PA)

%% System 2

clear all
% Load data
load S2_Single_EVPXI Action_Outcomes
load S2_Pairwise_EVPXI Action_DblOutcomes

% SINGLE PARAMETERS - EVPXI
Performance_single = squeeze((Action_Outcomes(2,:,:) - Action_Outcomes(1,:,:))./...
    (Action_Outcomes(3,:,:) - Action_Outcomes(1,:,:)))';
% Quantiles give us the range we display (max, min and average)
Q_single = quantile(Performance_single,[0.025 0.5 0.975]);

% COMBINATIONS OF PARAMETERS - EVPXI
Performance_double = squeeze((Action_DblOutcomes(2,:,:) - Action_DblOutcomes(1,:,:))./...
    (Action_DblOutcomes(3,:,:) - Action_DblOutcomes(1,:,:)))';
% Quantiles give us the range we display (max, min and average)
Q_double = quantile(Performance_double,[0.025 0.5 0.975]);

% Combine single parameter and combinations of parameters
PA = [Q_single Q_double];

% Write to csv
csvwrite('S2.csv', PA)

%% System 3

clear all
% Load data
load S3_Single_EVPXI Action_Outcomes
load S3_Pairwise_EVPXI Action_DblOutcomes

% SINGLE PARAMETERS - EVPXI
Performance_single = squeeze((Action_Outcomes(2,:,:) - Action_Outcomes(1,:,:))./...
    (Action_Outcomes(3,:,:) - Action_Outcomes(1,:,:)))';
% Quantiles give us the range we display (max, min and average)
Q_single = quantile(Performance_single,[0.025 0.5 0.975]);

% COMBINATIONS OF PARAMETERS - EVPXI
Performance_double = squeeze((Action_DblOutcomes(2,:,:) - Action_DblOutcomes(1,:,:))./...
    (Action_DblOutcomes(3,:,:) - Action_DblOutcomes(1,:,:)))';
% Quantiles give us the range we display (max, min and average)
Q_double = quantile(Performance_double,[0.025 0.5 0.975]);

% Combine single parameter and combinations of parameters
PA = [Q_single Q_double];

% Write to csv
csvwrite('S3.csv', PA)

%% System 4

clear all
load S4_Single_EVPXI Action_Outcomes
load S4_Pairwise_EVPXI Action_DblOutcomes

% SINGLE PARAMETERS - EVPXI
Performance_single = squeeze((Action_Outcomes(2,:,:) - Action_Outcomes(1,:,:))./...
    (Action_Outcomes(3,:,:) - Action_Outcomes(1,:,:)))';
% Quantiles give us the range we display (max, min and average)
Q_single = quantile(Performance_single,[0.025 0.5 0.975]);

% COMBINATIONS OF PARAMETERS - EVPXI
Performance_double = squeeze((Action_DblOutcomes(2,:,:) - Action_DblOutcomes(1,:,:))./...
    (Action_DblOutcomes(3,:,:) - Action_DblOutcomes(1,:,:)))';
% Quantiles give us the range we display (max, min and average)
Q_double = quantile(Performance_double,[0.025 0.5 0.975]);

% Combine single parameter and combinations of parameters
PA = [Q_single Q_double];

% Write to csv 
csvwrite('S4.csv', PA)
