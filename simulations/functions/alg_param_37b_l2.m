function np = alg_param_37b_l2(np)
% Step-size selection
% W. Ananduta
% 28/09/2020

np.A = cell(np.n,1);
a_add = 1;
% alphas
for i=1:np.n
    Ni=sum(np.Adj(i,:));
    
    % alphas
    
    
    a_pi = 1 + a_add;
    a_ch = 1 + a_add;
    a_ds = 1 + a_add;
    a_mg = 4 + np.n*np.d_mg + a_add; 
    
    a_tr =2.01*ones(Ni,1); 
    
    a_all = [a_pi;a_ch;a_ds;a_mg;a_tr;zeros(Ni,1)];
    a_ik = diag(a_all);
    np.A{i}= sparse(kron(eye(np.h),a_ik));    % Structure/order different than what's written on the paper, due to dynamics of storage


    

end


% beta
b = 0.49;
np.beta_tr = b*ones(np.n);

for y=1:np.b
    
    a_the = a_add;
    a_v = a_add;
    a_tg = 2 + a_add;
    a_p = (1 + a_add)*ones(length(np.B{y}),1);
    a_q = (0 + a_add)*ones(length(np.B{y}),1);
    
    a_all = [a_the; a_v; a_tg; a_p; a_q];
    np.A_DSO{y} = sparse(kron(diag(a_all),eye(np.h)));
    
    np.beta_pb(y) = 0.5/(length(np.B{y})+2*length(np.N_b{y})+a_add);
    

end

np.beta_tg = 1/(np.n+np.b+a_add);

np.gamma_mg = 1/(np.n+a_add);

end
