% Simulations
% extended P2P market, single simulation
% W. Ananduta
% 19/10/2020


clear all
close all
clc

ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = [1]; %uniform trading cost
%     
% run('case_37bus_N.m')
% % identify set of neighbors
% np.N = id_neigh(np.Adj);
% np.B = id_neigh(np.Adj_p);
% save(['case_jext_',date,'_B'],'np')

load(['case_jext_03-Feb-2021'])
%load( ['case_jext_',date,'_D1'])
% vary price of storage.
np.cstr = [0.02, 0.05, 0.5, 5]*100; % I_C
save(['case_jext_B_cst_12-Feb-2021'])
%% 
for cc = 1:length(np.cstr)
    
       
        
        np.st_un = ones(np.n,1);
        
        
        np = gen_param(np,ty);
        np = gen_cost(np,tc);
        
        np.q_st = np.cstr(cc)*ones(1,np.n);
        %np.c_st = np.cstr(cc)*ones(1,np.n);
        
        %save(['sim37b_B_sce_',num2str(sc),'_ct',num2str(ct(cc)),'_np'],'np')
        [s,sl,np] = sd_alg_eP2P_l1_f(np);
            
        % compute total cost
        [s,o] = com_cost(s,np);
        o.cst = np.c_st; 
       
%        save(['sim_jext_',date,'_l1_C1_',num2str(cc)],'o' )
        save(['sim_jext_B_cst_12-Feb-2021_',num2str(cc)],'o' )
        clearvars('s','sl')
%         %save(['sim37b_B_sce_',num2str(sc),'_ct',num2str(ct(cc)),'_np'],'np')
%         [s,sl,np] = sd_alg_eP2P_l2_np(np);
%             
%         % compute total cost
%         [s,o] = com_cost(s,np);
%         o.ct = ct;
%        
%         save(['sim_jext_',date,'_l2_C_',num2str(cc)],'o' )
        %clearvars('o','s')
        %if sc1 < 5
 %           save(['simJ1_n',num2str(sim.n_ag(sc)),'_tr_c',num2str(tr_cost(sc1))]);
        %else
        %    save(['sim4_n',num2str(sim.n_ag(sc)),'_tr_c_rand']);
        %end

end


%end
    

%end
% %%
% figure; plot(s.res);% set(gca,'yscale','log');
% %%
% figure;plot(s.ph_mg);