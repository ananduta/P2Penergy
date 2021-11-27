% Simulations
% extended P2P market, single simulation
% Varying coeff. of trading penalty
% W. Ananduta & G. Belgioioso
% 18/02/2021


clear all
close all
clc

% Add path of folder 'functions'
addpath([pwd,'\functions'])


ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = [1]; %uniform trading cost

% set the number of agents (EDIT)
n_agents = 10;
n_passive = 10;

% generate case
run('case_37bus_N.m')
% identify set of neighbors
np.N = id_neigh(np.Adj);
np.B = id_neigh(np.Adj_p);

% Selections of penalty coefficient
np.q_tr_var = [ 0.01, 0.03, 0.07, 0.1, 0.2, 1 ]*100;  % PARAMETERS VARIED

save(['case_sim_C_qtr_',date],'np')
%% 
for cc = 1:length(np.q_tr_var) 
    
    np.q_tr = np.qtr_var(cc);

    [s,sl,np] = sd_alg_eP2P_l1_f(np);

    % compute total cost
    [s,o] = com_cost(s,np);
    o.q_tr = np.q_tr; 
       
    save(['sim_C_qtr_',date,'_',num2str(cc)],'o')
    clearvars('s','sl');

end