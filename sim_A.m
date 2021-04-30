% Simulations
% extended P2P market, single simulation
% Line capacity violations
% W. Ananduta & G. Belgioioso
% 19/02/2020


clear all
close all
clc

% Add path of folder 'functions'
addpath([pwd,'\functions'])


ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = [1]; %uniform trading cost

% set the number of agents
n_agents = 120;
n_passive = 120;

% generate case
run('case_37bus_N.m')

% identify set of neighbors
np.N = id_neigh(np.Adj);
np.B = id_neigh(np.Adj_p);


% selections of line capacity constraints
sb = [((np.n+np.pas_ag)/np.b+2)*600 ((np.n+np.pas_ag)/np.b+2)*300]; % PARAMETERS VARIED

np.sb_set = sb;

save(['case_sim_A_',date],'np')
%% 
for cc = 1:length(sb)
        
    % set line capacity constraint
    np.s_bar = sb(cc)*ones(np.b);

    [s,sl,np] = sd_alg_eP2P_l1_f(np);

    % compute total cost
    [s,o] = com_cost(s,np);     

    save(['sim_A_',date,'_',num2str(cc)],'o' )
    clearvars('s','sl');

end