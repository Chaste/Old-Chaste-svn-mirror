function plot_s1s2_summary_graphs(data, All_APD90s, All_PCLs, MaxSlopes)
extra_concs_low = data.drug_low_ETPC;
extra_concs_high = data.drug_high_ETPC;

plot_against_DI = true; % As opposed to plotting against S2 which
                        % would also be (more!) sensible.

for i = 1:length(data.drug_indices_to_use)
    drug = data.drug_indices_to_use(i);
    pacing_cycle_lengths = All_PCLs{drug};
    APD90s_this_drug = All_APD90s{drug}; 
    MaxSlopes_this_drug = MaxSlopes{drug};
    
    concentration = data.drug_concentrations_text;
    if extra_concs_low(drug) > 0 % If data is missing it is given a '-1' in the file
        concentration{length(data.drug_concentrations_text)+1} = num2str(extra_concs_low(drug));
        concentration{length(data.drug_concentrations_text)+2} = num2str((extra_concs_low(drug)+extra_concs_high(drug))/2.0);
        concentration{length(data.drug_concentrations_text)+3} = num2str(extra_concs_high(drug));
    end

    num_cols_on_big_fig = 4;
    fignum=ceil(drug/num_cols_on_big_fig);

    % Control plot

    for conc = 1:length(concentration)
        pacing_rates=[];
        APDs=[];
        stan_sol=[];
        size(APD90s_this_drug)
        APD90s = APD90s_this_drug{conc}';
        temp = pacing_cycle_lengths{conc};
        s2s = temp(:,1);
        DIs = temp(:,2);
        if conc==1
            if plot_against_DI
                control = [APD90s DIs];
            else
                control = [APD90s s2s];
            end
        else
            figure(200+fignum)
            subplot(10,num_cols_on_big_fig,(conc-2)*num_cols_on_big_fig + drug - (fignum-1)*num_cols_on_big_fig)
            plot(control(:,2),control(:,1),'k.-')
            hold on
            if plot_against_DI
                plot(DIs,APD90s,'r.-')
            else
                plot(s2s,APD90s,'r.-')
            end
            title([data.drug_names{drug} ' at ' concentration{conc} ' nM, Max Slope = ' num2str(MaxSlopes_this_drug(conc))])
            xlim([0 1000])
            if max(APD90s)>400
                ylim([0 max(APD90s)])
            else
                ylim([0 400])
            end
        end
    end
end