function [max_slope alternans_start separate_apd_frequency areas APD90s pacing_cycle_lengths] = LoadDynamicRestitutionData( drug, data )
%LoadDynamicRestitutionData Load the dynamic protocol results
%   We find the APDs and detect alternans onset.

% Dynamic protocol details from Chaste Simulation
num_pulses_to_calc_APD = 8;
pacing_cycle_lengths = [1000:-100:600 550 500 480:-20:360 350:-10:300 295:-5:100]';

% Does this drug have clinical dose information 
% and therefore extra conc sims?
extra_concs_low = data.drug_low_ETPC;
extra_concs_high = data.drug_high_ETPC;
if extra_concs_low(drug) > 0 % If data is missing it is given a '-1' in the file
    data.drug_concentrations_text{9} = num2str(extra_concs_low(drug));
    data.drug_concentrations_text{10} = num2str((extra_concs_low(drug)+extra_concs_high(drug))/2.0);
    data.drug_concentrations_text{11} = num2str(extra_concs_high(drug));
    data.drug_concentrations_text{12} = num2str(10.0*extra_concs_high(drug));
end

filename = [data.foldername 'Dynamic/voltage_results.dat'];
all_data = importfile(filename);
dynamic_data = all_data.data;
stan_sol=[];

% Set up correct sized array
alternans_control_idx = length(pacing_cycle_lengths);
APD90s{length(data.drug_concentrations_text),length(pacing_cycle_lengths),num_pulses_to_calc_APD} = [];

for conc = 1:length(data.drug_concentrations_text)
    for pacing_index = 1:length(pacing_cycle_lengths)
        for i=1:length(dynamic_data)
            % Skim through the data file...
            if (dynamic_data(i,1) == drug-1)
                % If this is our drug of interest
                if (abs(dynamic_data(i,2)-str2double(data.drug_concentrations_text{conc}))<1e-6)
                    % If this is the concentration of interest
                    if (abs(dynamic_data(i,3)-pacing_cycle_lengths(pacing_index))<1e-6)
                        % If this is our PCL of interest
                        for beat=1:num_pulses_to_calc_APD
                            APD90s{conc, pacing_index, beat} = dynamic_data(i, 3+beat);
                        end
                    end
                end
            end
        end        
    end
    disp(['Conc = ' num2str(data.drug_concentrations_text{conc}) ' nM'])
    
    pacing_cell = [];
    apd_cell = [];
    % Assemble vectors for each beat
    for beat_index=1:num_pulses_to_calc_APD
        pacing_rates=[];
        APDs=[];
        for pacing_index = 1:length(pacing_cycle_lengths)
            if APD90s{conc, pacing_index, beat_index}>0
                pacing_rates = [pacing_rates; pacing_cycle_lengths(pacing_index)];
                APDs = [APDs; APD90s{conc, pacing_index, beat_index}];
                % Record a standard (drug free) solution so we can plot
                % it as a control...
                if (data.drug_concentrations_text{conc}=='0')
                    stan_sol = [stan_sol; pacing_cycle_lengths(pacing_index) APD90s{1, pacing_index, beat_index} ];
                end
            end
        end
        if (conc==1)
            standard_solution{beat_index} = stan_sol;
            stan_sol=[];
        end
        
        figure(8523);
        plot(pacing_rates,APDs,'b.-')
        if max(APDs)>400
            ylim([0 max(APDs)])
        else
            ylim([0 400])
        end
        xlim([min(pacing_cycle_lengths) max(pacing_cycle_lengths)]);
        hold on
        pacing_cell{beat_index} = pacing_rates;
        apd_cell{beat_index} = APDs;
    end
        
    [alternans_start{conc} separate_apd_frequency{conc} areas(conc) alternans_start_idx max_slope(conc)]= detect_alternans_starts(pacing_cell,apd_cell,num_pulses_to_calc_APD,pacing_cycle_lengths,standard_solution, alternans_control_idx);
    
    if conc==1
        alternans_control_idx=alternans_start_idx;
    end
       
    figure(8523)
    xlabel('Pacing Cycle Length (ms)')
    xlim([min(pacing_cycle_lengths) max(pacing_cycle_lengths)]);
    ylabel('APD90')
    title([data.model_name ' ' data.drug_names{drug} ' at ' ...
        data.drug_concentrations_text{conc} 'nM, Dynamic Restitution' ...
        ' Max slope is ' num2str(max_slope(conc)) ])
    hold on
end
hold off

function [alternans_start_freq ead_onset_freq area_between alternans_start_idx max_slope] = detect_alternans_starts(pacing_rates, apds, num_pulses_taken, pacing_cycle_lengths,standard_solution,alternans_control_idx)
alternans_start_freq = -1;
alternans_start_idx = length(pacing_cycle_lengths);
ead_onset_freq = -1;

% Initially figure out where there are still "num_pulses_taken" separate
% APDs
for pacing_rate=length(pacing_cycle_lengths):-1:1
    num_at_this_freq = 0;
    for beat_index=1:num_pulses_taken
        pacing_vec = pacing_rates{beat_index};
        apd_vec = apds{beat_index};
        if (find(abs(pacing_vec-pacing_cycle_lengths(pacing_rate))<1e-6))
            num_at_this_freq = num_at_this_freq+1;
        end
    end
    if (num_at_this_freq==num_pulses_taken)
        if pacing_rate~=length(pacing_cycle_lengths)
            ead_onset_freq = pacing_cycle_lengths(pacing_rate+1);
        end
        all_separate_apds_freq_idx = pacing_rate;
        break;
    end
end


% Now look for alternans onset
% Work out mean deviation before onset and look for
% a big increase...
standard_deviation=[];
for pacing_rate=1:all_separate_apds_freq_idx
    num_at_this_freq = 0;
    apds_this_freq = [];
    for beat_index=1:num_pulses_taken
        pacing_vec = pacing_rates{beat_index};
        apd_vec = apds{beat_index};
        apds_this_freq(beat_index) = apd_vec(pacing_rate);
    end
    standard_deviation(pacing_rate) = std(apds_this_freq);
    
    if standard_deviation(end)>1 % We want the standard deviation to be over 1 ms.
        alternans_start_freq = pacing_cycle_lengths(pacing_rate);
        alternans_start_idx = pacing_rate;
        break;
    end
end

% Work out mean control and drug lines before alternans
for beat_idx = 1:num_pulses_taken
    temp = standard_solution{beat_idx};
    temp2 = apds{beat_idx};
    if alternans_start_idx < alternans_control_idx
       sum_to_here =  alternans_start_idx-1;
    elseif all_separate_apds_freq_idx < alternans_control_idx
       sum_to_here = all_separate_apds_freq_idx;
    else
       sum_to_here =  alternans_control_idx-1;
    end
    control_vec(:,beat_idx) = temp(1:sum_to_here,2);
    drug_vec(:,beat_idx) = temp2(1:sum_to_here);
end

assert(length(control_vec)==length(drug_vec));
control_vec = mean(control_vec,2);
drug_vec = mean(drug_vec,2);

% We can only provide most of the outputs if there is a regime with no 
% alternans. If alternans occur from the start (low freq pacing) then
% return zeros for these measures.
if length(control_vec)>1
    slopes = diff(drug_vec(1:sum_to_here))./diff(pacing_cycle_lengths(1:sum_to_here));
    [val, idx] = max(abs(slopes));
    max_slope = slopes(idx);
    % Area between dynamic pacing under drug and control.
    area_between = trapz(pacing_cycle_lengths(1:length(control_vec)),control_vec-drug_vec);
else
    max_slope = 0;
    area_between = 0;
end


