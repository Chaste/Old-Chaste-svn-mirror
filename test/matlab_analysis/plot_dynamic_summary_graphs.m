function plot_dynamic_summary_graphs(data, APD90s, PCLs)
extra_concs_low = data.drug_low_ETPC;
extra_concs_high = data.drug_high_ETPC;

for i = 1:length(data.drug_indices_to_use)
    drug = data.drug_indices_to_use(i);
    pacing_cycle_lengths = PCLs{drug}; % (Same for every drug??)
    APD90s_this_drug = APD90s{drug}; % (Same for every drug??)
    
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
        for beat_index=1:4
            pacing_rates=[];
            APDs=[];
            stan_sol=[];
            for pacing_index = 1:length(pacing_cycle_lengths)
                if APD90s_this_drug{conc, pacing_index, beat_index}>0
                    pacing_rates = [pacing_rates; pacing_cycle_lengths(pacing_index)];
                    APDs = [APDs; APD90s_this_drug{conc, pacing_index, beat_index}];
                    % Record a standard (drug free) solution so we can plot
                    % it as a control...
                    if conc==1
                        stan_sol = [stan_sol; pacing_cycle_lengths(pacing_index) APD90s_this_drug{1, pacing_index, beat_index} ];
                        standard_solution{beat_index} = stan_sol;
                    end
                end
            end
            if conc>1
                figure(300+fignum)
                subplot(10,num_cols_on_big_fig,(conc-2)*num_cols_on_big_fig + drug - (fignum-1)*num_cols_on_big_fig)
                stan_sol = standard_solution{beat_index};
                plot(stan_sol(:,1),stan_sol(:,2),'k.-')
                hold on
                plot(pacing_rates,APDs,'r.-')
                title([data.drug_names{drug} ' at ' concentration{conc} ' nM'])
                xlim([min(pacing_cycle_lengths) max(pacing_cycle_lengths)]);
                if max(APDs)>400
                    ylim([0 max(APDs)])
                else
                    ylim([0 400])
                end
            end
        end
    end
end
