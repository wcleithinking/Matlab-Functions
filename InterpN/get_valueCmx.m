function Cmx = get_valueCmx(Cmx_row,...
    Ma_list,alpha_deg_list,beta_deg_list, ...
    deltaphi_deg_list,deltapsi_deg_list,deltagamma_deg_list, N, ...
    Ma,alpha_deg,beta_deg, ...
    deltaphi_deg,deltapsi_deg,deltagamma_deg)

%   Get the value of aerodynamic coefficient Cmx based on the following inputs:
%       1) Cmx_row;
%       2) Ma_list,alpha_deg_list,beta_deg_list,deltaphi_deg_list,detapsi_deg_list,deltagamma_deg_list,N;
%       3) Ma,alpha_deg,beta_deg,deltaphi_deg,deltapsi_deg,deltagamma_deg.

%% 1. Find the left indexes of inputs
Ma_leftindex = get_leftindex(Ma_list,Ma);
alpha_deg_leftindex = get_leftindex(alpha_deg_list,alpha_deg);
beta_deg_leftindex = get_leftindex(beta_deg_list,beta_deg);
deltaphi_deg_leftindex = get_leftindex(deltaphi_deg_list,deltaphi_deg);
deltapsi_deg_leftindex = get_leftindex(deltapsi_deg_list,deltapsi_deg);
deltagamma_deg_leftindex = get_leftindex(deltagamma_deg_list,deltagamma_deg);

%% 2. get the value of Cmx(Ma,alpha,beta,deltaphi,deltapsi,deltagamma)
Cmx = Cmx_row( prod(N(2:6))*(Ma_leftindex-1) + prod(N(3:6))*(alpha_deg_leftindex-1) + ...
    prod(N(4:6))*(beta_deg_leftindex-1) + prod(N(5:6))*(deltaphi_deg_leftindex-1) + ...
    prod(N(6))*(deltapsi_deg_leftindex-1) + deltagamma_deg_leftindex);
