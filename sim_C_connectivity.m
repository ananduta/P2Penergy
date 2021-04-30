% Simulations
% extended P2P market, single simulation
% Varying connectivity in the trading network
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

% set connectivity
np.con = [ 0.2:0.2:1 ]; % PARAMETERS VARIED

save(['case_sim_C_con_',date],'np')
%% 
for cc = 1:length(np.con)
        
        % generate adjacency matrix
        np.Adj = sparse(randconG(np.n,np.con(cc)));
        
        % update parameters of the cost and constraints
        np = gen_param(np,ty); 
        np = gen_cost(np,tc); 
        np.N = id_neigh(np.Adj);
        
        %
        [s,sl,np] = sd_alg_eP2P_l1_f(np);
            
        % compute total cost
        [s,o] = com_cost(s,np);
        o.con = np.con(cc);
        o.Adj = np.Adj;

        save(['sim_jext_C_con_',date,'_',num2str(cc)],'o' )
        clearvars('s','sl')

end