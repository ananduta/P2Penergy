% Generate a case using IEEE 14-bus network

% set time horizon
np.h= 24;


%% Adjacency matrix of trading network
np.n = 10;
np.Adj = sparse(randconG(np.n,0.6));

%% Adjacency matrix of physical network
np.b = 36;
np.Adj_p = zeros(np.b);
np.Bnet = zeros(np.b);
np.Gnet = zeros(np.b);

% operating voltage
np.v_op = 4.8^2*1000; %kV

load('case_37bus.mat');
tab = case_37bus;
for i= 1:length(tab(:,1))
    np.Adj_p(tab(i,1),tab(i,2)) = 1;
    
    if  tab(i,4) == 721
        z_b = 0.2926 + 0.1973i;
    elseif tab(i,4) == 722
        z_b = 0.4751 + 0.2973i;
    elseif tab(i,4) == 723
        z_b = 1.2936 + 0.6713i;
    elseif tab(i,4) == 724
        z_b = 2.0952 + 0.7758i;
    end
    z= z_b*tab(i,3)/5280;
    np.Bnet(tab(i,1),tab(i,2)) = abs(imag(1/z));
    np.Gnet(tab(i,1),tab(i,2)) = abs(real(1/z));
end
np.Adj_p = sparse(np.Adj_p + np.Adj_p');
np.Bnet = sparse(np.Bnet + np.Bnet');
np.Bnet = np.v_op*np.Bnet;
np.Gnet = sparse(np.Gnet + np.Gnet');
np.Gnet = np.v_op*np.Gnet;

% index of bus where each prosumer is connected to
%np.B_n = idx_bus(np.N_b,np.n);
np.B_n = randperm(np.b,mod(np.n,np.b));
for i = 1:floor(np.n/np.b)
    np.B_n = [np.B_n randperm(np.b,np.b)];
end

%np.B_n = [randperm(np.b,np.b) randperm(np.b,np.b) randperm(np.b,100-72)];

% set of prosumers attached to the bus
for y=1:np.b
    np.N_b{y}= find(np.B_n==y);
end

% Busses connected to the main grid
np.B_mg = [36];

%% assign components to each node
n = np.n;
b = np.b;
% type of load
% randomly assign the type of load profile   
np.t_lpr = randi([1 6], n+b,1);
% 0 = no load
%np.t_lpr(n+1) = 0;
%np.t_lpr(n+7) = 0;
%np.t_lpr(n+8) = 0;
    
% storage units (no storage units in [Sousa, et. al.,2019] and [Le Cadre, et al, 2019]
np.st_un = zeros(n,1);
np.st_un = randi([0 1], n,1);

% dispatchable units
%np.d_un = randi([0 1], n,1);
%np.d_un = [ones(1, 10) zeros(1,n-10)];

 n_di=4;
 ag_dg = randperm(np.n,n_di); 
 np.d_un = zeros(n,1);
 for i = 1:n_di
     np.d_un(ag_dg(i)) = 1;
end
%np.d_un = zeros(n,1);

% PV generation units
np.r_un = zeros(n+b,1);
for i=1:n
    if np.st_un == 1
        np.r_un(i) = randi([0 4]);
    end
end
%np.r_un = [randi([1 4], n,1); zeros(b,1)];


% assign parameters in the local constraints
np = gen_param(np,ty); 

% assign per-unit costs
np = gen_cost(np,tc); %MODIFY

% generate load and non-dispatchable profiles
[np.Pd,np.Pl,np.Pr] = gen_Pd_ext(np.n,np.t_lpr,np.r_un,np.b);
%[np.Pd,np.Pl,np.Pr] = gen_Pd_ext1(np.n,np.t_lpr,np.r_un,np.b);
%np.h = size(np.Pd,2); % time horizon

% initial condition of variables
np.init = 0;
