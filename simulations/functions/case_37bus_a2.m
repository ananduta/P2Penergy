% Generate a case using IEEE 14-bus network

% load base case (Power demand and Adjacency matrix)
load('37b_base.mat')

% set time horizon
np.h= 24;


%% Adjacency matrix of trading network
np.n = 25+36;
np.Adj = min([Adj zeros(25,36); zeros(36,25+36)]+randconG(np.n,0.6),1);

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
np.Adj_p = np.Adj_p + np.Adj_p';
np.Bnet = np.Bnet + np.Bnet';
np.Bnet = np.v_op*np.Bnet;
np.Gnet = np.Gnet + np.Gnet';
np.Gnet = np.v_op*np.Gnet;

% set of prosumers attached to the bus
np.N_b{1} = 7;
np.N_b{2} = 1;
np.N_b{3} = [2,3];
np.N_b{4} = double.empty;
np.N_b{5} = double.empty;
np.N_b{6} = 8;
np.N_b{7} = 4;
np.N_b{8} = [5,6];
np.N_b{9} = double.empty;
np.N_b{10} = double.empty;
np.N_b{11} = double.empty;
np.N_b{12} = double.empty;
np.N_b{13} = 9;
np.N_b{14} = double.empty;
np.N_b{15} = double.empty;
np.N_b{16} = double.empty;
np.N_b{17} = double.empty;
np.N_b{18} = 10;
np.N_b{19} = double.empty;
np.N_b{20} = double.empty;
np.N_b{21} = 11;
np.N_b{22} = 12;
np.N_b{23} = [13,14];
np.N_b{24} = double.empty;
np.N_b{25} = double.empty;
np.N_b{26} = 15;
np.N_b{27} = 16;
np.N_b{28} = [17,18];
np.N_b{29} = double.empty;
np.N_b{30} = 19;
np.N_b{31} = 20;
np.N_b{32} = 21;
np.N_b{33} = [22,23];
np.N_b{34} = double.empty;
np.N_b{35} = [24,25];
np.N_b{36} = double.empty;

for y = 1:np.b
    np.N_b{y} = [np.N_b{y}, 25+y];
end
% index of bus where each prosumer is connected to
np.B_n = idx_bus(np.N_b,np.n);

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

% dispatchable units
%np.d_un = randi([0 1], n,1);
np.d_un = [d_un; zeros(b,1)];

% PV generation units
%np.r_un = [randi([1 4], n,1); zeros(b,1)];


% assign parameters in the local constraints
np = gen_param(np,ty); 

% assign per-unit costs
np = gen_cost(np,tc); %MODIFY

% generate load and non-dispatchable profiles
%[np.Pd,np.Pl,np.Pr] = gen_Pd_ext(np.n,np.t_lpr,np.r_un,np.b);
%np.h = size(np.Pd,2); % time horizon
np.Pd = [Pd; zeros(np.b,np.h)];
