% Simulations
% extended P2P market
% Varying connectivity (vs. #iterations and others)
% W. Ananduta & G. Belgioioso
% 18/02/2021


%clear all
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
run('case_37bus_N.m')
% identify set of neighbors
np.N = id_neigh(np.Adj);
np.B = id_neigh(np.Adj_p);

% set connectivity
np.con = [ 0.1:0.1:1 ];% DEFINE SPARSITY OF THE TRADING NETWORK 
n_sim =10; % DEFINE THE NUMBER OF SIMULATIONS PER CASE
%% 
for cc = 10:length(np.con)
    for cc1 = 2:5
        %% Adjacency matrix of trading network
        np.Adj = randconG(np.n,np.con(cc));
        
        np.N = id_neigh(np.Adj);

        %semi-decentralized
        [s,sl,np] = sd_alg_eP2P_l1_f(np);

        % compute total cost
        [s,o] = com_cost(s,np);
        
        r.er_max = np.er_max;
        r.iter_tab(cc,cc1) = o.iter; 
        r.J_tab(cc,cc1) = o.Jt;
        r.Jal_tab(cc,cc1) = o.Jall;
        
        % compute total trade
        Ptr_t = zeros(np.h,1);
    
        for i=1:np.n
            for j = 1:np.n
                if ~isempty(o.p_tr{i,j}) 
                    for k=1:np.h
                        Ptr_t(k,1) = Ptr_t(k,1) + max(0,o.p_tr{i,j}(k,end));
                    end
                end
            end
        end
        r.Ptr_t(cc,cc1) = sum(Ptr_t);
        
        % save results
        save(['sim_D_con_',date],'r')
        clearvars('s','sl')
    end

end