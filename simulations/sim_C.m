% Simulations
% extended P2P market, single simulation
% Varying cost coeff. of trading + disallowed trading
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

% selections of cost coefficients.
ctr_var = [ 0.02, 0.05, 0.07, 0.08, 0.12, 0.2, 0.5, 1 ]*100; 
np.ctr_var = ctr_var;

save(['case_sim_C_',date],'np')
%% 
for cc = 1:length(ctr_var)+1
    
    np.c_tr = np.ctrl(cc)*np.Adj;
    
    if cc == length(ctr_var)+1
        np.c_tr = np.ctrl(3)*np.Adj;
        np.pt_max = 0*np.Adj;
    end
    
    %
    [s,sl,np] = sd_alg_eP2P_l1_f(np);

    % compute total cost
    [s,o] = com_cost(s,np);
    
    %
    o.c_tr = np.c_tr; 
    o.pt_max = np.pt_max;

    save(['sim_C_',date,'_',num2str(cc)],'o')
    clearvars('s','sl');
       
end


%end
    

%end
% %%
% figure; plot(s.res);% set(gca,'yscale','log');
% %%
% figure;plot(s.ph_mg);