clc; clear; close all;
Ma_list = 1.5:0.5:6;
alpha_deg_list = -25:5:25;
beta_deg_list = -25:5:25;
deltaphi_deg_list = -30:10:30;
deltapsi_deg_list = -30:10:30;
deltagamma_deg_list = -30:10:30;
deg2rad = pi/180;
alpha_rad_list = alpha_deg_list*deg2rad;
beta_rad_list  = beta_deg_list*deg2rad;
deltaphi_rad_list = deltaphi_deg_list*deg2rad;
deltapsi_rad_list = deltapsi_deg_list*deg2rad;
deltagamma_rad_list = deltagamma_deg_list*deg2rad;
N_Ma = length(Ma_list);
N_alpha = length(alpha_rad_list);
N_beta = length(beta_rad_list);
N_deltaphi = length(deltaphi_rad_list);
N_deltapsi = length(deltapsi_rad_list);
N_deltagamma = length(deltagamma_rad_list);
N = [N_Ma,N_alpha,N_beta,N_deltaphi,N_deltapsi,N_deltagamma];
p = [1e-3,0.2,-0.2, ...
    0.15,-0.5,0.5, ...
    0.24,0.05,0.05];
Cmx_row = zeros(prod(N),1);
counter = 1;
for i = 1:N(1)
    for j = 1:N(2)
        for k = 1:N(3)
            for l = 1:N(4)
                for m = 1:N(5)
                    for n = 1:N(6)
                        Cmx_row(counter) = p(1)+p(2)*alpha_rad_list(j)+p(3)*beta_rad_list(k)+ ...
                            p(4)*alpha_rad_list(j)*beta_rad_list(k)+ ...
                            p(5)*beta_rad_list(k)*deltaphi_rad_list(l)+ ...
                            p(6)*alpha_rad_list(j)*deltapsi_rad_list(m)+...
                            ( p(7)+p(8)*alpha_rad_list(j)^2+p(9)*beta_rad_list(k)^2)*deltagamma_rad_list(n);
                        counter = counter + 1;
                    end
                end
            end
        end
    end
end
save data_Cmx.mat Cmx_row Ma_list alpha_deg_list beta_deg_list deltaphi_deg_list deltapsi_deg_list deltagamma_deg_list N