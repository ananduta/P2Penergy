% Simulations
% extended P2P market, single simulation
% Peak-shaving effect of storage units
% W. Ananduta & G. Belgioioso
% 18/02/2021


clear all
close all
clc

% Add path of folder 'functions'
addpath([pwd,'\functions'])


ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = [1]; %uniform trading cost

% set the number of agents
n_agents = 10;
n_passive = 10;

% generate case
run('case_37bus_N_B.m')
% identify set of neighbors
np.N = id_neigh(np.Adj);
np.B = id_neigh(np.Adj_p);

save(['case_sim_B_',date],'np')
%% 
for cc = 1:3
    
    if cc == 2
        % half of agents have storage units
        n_st=floor(np.n/2);
         ag_st = randperm(np.n,n_st); 
         np.st_un = zeros(n,1);
         for i = 1:n_st
             np.st_un(ag_st(i)) = 1;
         end
    elseif cc == 3
        % all agents have storage units
        np.st_un = ones(np.n,1);
    end

    % update constraints and cost parameters
    np = gen_param(np,ty);
    np = gen_cost(np,tc);
    
    % run distributed algorithm
    [s,sl,np] = sd_alg_eP2P_l1_f(np);

    % compute total cost
    [s,o] = com_cost(s,np);     

    % record possesion of storage unit of each agent
    o.st_un = np.st_un;

    save(['sim_B_',date,'_',num2str(cc)],'o')
    clearvars('s','sl');


end