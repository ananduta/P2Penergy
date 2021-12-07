function np = gen_param(np,ty,Ts)
% Set the parameters in the constraints
% W. Ananduta
% 22/07/2019
% Edit for extension: 30/09/2020

% Inputs:
% np = struct
% ty = type of case

% Output:
% np = struct

if nargin < 3
    Ts=1;
end

% homogenous case study, ty==1
% local constraints
np.pdg_min = zeros(1,np.n);
np.pdg_max = 20*ones(1,np.n);
np.p_ds = 5*ones(1,np.n);
np.p_ch = 5*ones(1,np.n);
np.x_min = 20*ones(1,np.n);
np.x_max = 100*ones(1,np.n);

% storage unit
T = 60*60*Ts; % s
J2Wh = 1/3600;
np.a_st = ones(1,np.n);
np.b_st = T*J2Wh*ones(1,np.n);
np.x0 = np.x_min+30;

np.eta_ch = 0.95*ones(1,np.n);
np.eta_ds = 0.94*ones(1,np.n);

% power transfer constraints
np.pt_max = 30*np.Adj;

% power to main grid constraints
np.pmg_min = 0;
np.pmg_max = (np.n+np.b)*300;

np.pmg_min_a = 0;
%np.pmg_max_a = 300;

% heterogenous case study, ty==0
if ty==0
    for i=1:np.n
        if np.d_un(i) ==1
            np.pdg_max(i) = np.t_lpr(i)^2*20;
        else
            np.pdg_max(i) = 0;
        end
        if np.st_un(i) == 0
            np.p_dh(i) = 0;
            np.p_ch(i) = 0;
            np.x0(i) = 0;
            np.x_min(i) = 0;
            np.x_max(i) = 0;
        else
            np.p_dh(i) = 10*np.t_lpr(i);
            np.p_ch(i) = 10*np.t_lpr(i);
            np.x_max(i) =100*(1+np.t_lpr(i)/max(np.t_lpr));
        end
            
    end
end

% parameters of the physical network
np.theta_max = 15/180*pi*ones(np.b,1);
np.theta_min = -15/180*pi*ones(np.b,1);

np.v_max = 1.05*ones(np.b,1);
np.v_min = 0.95*ones(np.b,1);

for yy = 1:length(np.B_mg)
    y = np.B_mg(yy);
    np.theta_max(y) = 0;
    np.theta_min(y) = 0;
    np.v_max(y) = 1;
    np.v_min(y) = 1;
end
np.s_bar = (np.n+np.pas_ag)*100*ones(np.b);
%np.s_bar = 80*ones(np.b);

end