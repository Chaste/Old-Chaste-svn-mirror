function [max_slope_this_drug APD90s_all_concs s2s_all_concs] = LoadS1S2data( drug, input_data )
%LoadS1S2data Load the S1-S2 protocol results
%   We analyse the results and return APD90s and Diastolic intervals.

% Does this drug have clinical dose information 
% and therefore extra conc sims?
extra_concs_low = input_data.drug_low_ETPC;
extra_concs_high = input_data.drug_high_ETPC;
if extra_concs_low(drug) > 0 % If data is missing it is given a '-1' in the file
    input_data.drug_concentrations_text{9} = num2str(extra_concs_low(drug));
    input_data.drug_concentrations_text{10} = num2str((extra_concs_low(drug)+extra_concs_high(drug))/2.0);
    input_data.drug_concentrations_text{11} = num2str(extra_concs_high(drug));
    input_data.drug_concentrations_text{12} = num2str(10.0*extra_concs_high(drug));
end

% Load the steady-state zero drug case to compare and plot with.
file_name = [input_data.foldername 'S1S2/voltage_results.dat']

all_data = importfile(file_name);
headers = all_data.colheaders;
data = all_data.data;

% We are interested in this drug...
for conc = 1:length(input_data.drug_concentrations_text)
    
    clear APD90s
    clear s2s
    clear DIs
    
    s2_counter=1;
    min_s2 = 1e10; % Check we don't count an S2 twice if it is used as a dose too.
    for i=1:length(data)
        % Skim through the data file...
        
        if data(i,1) == drug-1
            % If this is our drug of interest
            
            if (abs(data(i,2)-str2double(input_data.drug_concentrations_text{conc}))<1e-6)
                % If this is the concentration of interest
                
                if data(i,3) < min_s2 
                    % (Only record monotonically decreasing S2s to avoid 
                    % duplicating data where a standard conc = a dose)
                    s2s(s2_counter) = data(i,3);
                    LastButOneAPD90 = data(i,5);
                    APD90s(s2_counter) = data(i,7);

                    s2_counter = s2_counter+1;
                    min_s2 = data(i,3);
                end
            end
        end        
    end
    DIs = s2s-LastButOneAPD90;
    slopes = diff(APD90s)./diff(DIs);
    [~, idx] = max(abs(slopes));
    max_slope_this_drug(conc) = slopes(idx);
    
    figure
    subplot(2,1,1)
    plot(s2s,APD90s,'b.-')
    xlabel('S2 period')
    ylabel('APD90')
    title([input_data.model_name ' ' input_data.drug_names{drug} ' at ' input_data.drug_concentrations_text{conc} 'nM'])
    
    subplot(2,1,2)
    plot(DIs,APD90s,'b.-')
    title(['S1S2 Max slope is ' num2str(max_slope_this_drug(conc)) ])
    ylabel('APD90')
    xlabel('Diastolic Interval')
    
    APD90s_all_concs{conc} = APD90s;
    s2s_all_concs{conc} = [s2s' DIs'];
end




