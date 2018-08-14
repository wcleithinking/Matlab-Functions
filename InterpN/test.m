clc;clear; close all;
generate_data;
counter = 1;
Cmx_fit = zeros(prod(N),1);
error = zeros(prod(N),1);
for i = 1:N(1)
    for j = 1:N(2)
        for k = 1:N(3)
            for l = 1:N(4)
                for m = 1:N(5)
                    for n = 1:N(6)
                        Cmx_fit(counter) = fit_valueCmx(Cmx_row, ...
                            Ma_list,alpha_deg_list,beta_deg_list, ...
                            deltaphi_deg_list,deltapsi_deg_list,deltagamma_deg_list, N, ...
                            Ma_list(i),alpha_deg_list(j),beta_deg_list(k),...
                            deltaphi_deg_list(l),deltapsi_deg_list(m),deltagamma_deg_list(n));
                        error(counter) = Cmx_fit(counter)-Cmx_row(counter);
                        counter = counter + 1;
                    end
                end
            end
        end
    end
end
plot(error);