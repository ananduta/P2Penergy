% Simulations
% extended P2P market
% Varying number of agents
% W. Ananduta & G. Belgioioso
% 18/02/2021


%clear all
close all
clc

% Add path of folder 'functions'
addpath([pwd,'\functions'])


ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = [1]; %uniform trading cost
%     

n_ag =[40:10:90]; % DEFINE THE NUMBER OF AGENTS HERE AS AN ARRAY
sp = 0.6; % DEFINE SPARSITY OF THE TRADING NETWORK 
n_sim =10; % DEFINE THE NUMBER OF SIMULATIONS PER CASE
%% 
for cc = 6:length(n_ag)
    for cc1 = 1:6
        
        n_agents = n_ag(cc);
        n_passive = 80;
        con = sp;
        
        run('case_37bus_N_D.m')
        np.N = id_neigh(np.Adj);
        np.B = id_neigh(np.Adj_p); 

        %semi-decentralized
        [s,sl,np] = sd_alg_eP2P_l1_f(np);

        % compute total cost
        [s,o] = com_cost(s,np);

        r.er_max = np.er_max;
        r.iter_tab(cc,cc1) = o.iter; 
        r.J_tab(cc,cc1) = o.Jt;
        r.Jal_tab(cc,cc1) = o.Jall;
        
        save(['sim_D_ag_',date],'r')
        clearvars('s','sl')
        
    end

end

for cc = 6:length(n_ag)
    for cc1 = 1:n_sim
        
        n_agents = n_ag(cc);
        n_passive = 80;
        con = sp;
        
        run('case_37bus_N_D.m')
        np.N = id_neigh(np.Adj);
        np.B = id_neigh(np.Adj_p); 

        %semi-decentralized
        [s,sl,np] = sd_alg_eP2P_l1_f(np);

        % compute total cost
        [s,o] = com_cost(s,np);

        r.er_max = np.er_max;
        r.iter_tab(cc,cc1) = o.iter; 
        r.J_tab(cc,cc1) = o.Jt;
        r.Jal_tab(cc,cc1) = o.Jall;
        
        save(['sim_D_ag_',date],'r')
        clearvars('s','sl')
        
    end

end
    