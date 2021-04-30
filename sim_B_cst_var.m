% Simulations
% extended P2P market, single simulation
% Varying cost coeff. of using storage units
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
n_agents = 50;
n_passive = 80;

% generate case
run('case_37bus_N_B_cst_var.m')
% identify set of neighbors
np.N = id_neigh(np.Adj);
np.B = id_neigh(np.Adj_p);

% selections of cost coefficients.
cst_var = [0.01, 0.05, 0.1, 0.2 ]*100; % EDIT
np.cst_var = cst_var;

save(['case_sim_B_',date],'np')
%% 
for cc = 1:length(cst_var)
    
    
    % update constraints and cost parameters
    np = gen_param(np,ty);
    
    np.q_st = cst_var(cc)*ones(1,np.n);
    
    % run distributed algorithm
    [s,sl,np] = sd_alg_eP2P_l1_f(np);

    % compute total cost
    [s,o] = com_cost(s,np);     

    % record possesion of storage unit of each agent
    o.q_st = np.q_st;

    save(['sim_B_cst_',date,'_',num2str(cc)],'o')
    clearvars('s','sl');


end