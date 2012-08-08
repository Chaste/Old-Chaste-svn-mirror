function [AP_triangulation APD90s APD50s Vmax Upstroke] = LoadSteadyStateData(drug, input_data, filename)

input_data.foldername
% Load the steady-state zero drug case to compare and plot with.
file_name = [input_data.foldername ...
    num2str(input_data.hertz) 'Hz_SteadyState/'...
    filename]

all_data = importfile(file_name);
headers = all_data.colheaders;
data = all_data.data;

drug_c_index = drug -1;
start_idx = min(find(abs(data(:,1)-drug_c_index)<1e-6));

drug_concentrations = input_data.drug_concentrations_text;
extra_concs_low = input_data.drug_low_ETPC;
extra_concs_high = input_data.drug_high_ETPC;
if extra_concs_low(drug) > 0 % If data is missing it is given a '-1' in the file
    drug_concentrations{9} = num2str(extra_concs_low(drug));
    drug_concentrations{10} = num2str((extra_concs_low(drug)+extra_concs_high(drug))/2.0);
    drug_concentrations{11} = num2str(extra_concs_high(drug));
    drug_concentrations{12} = num2str(10*extra_concs_high(drug));
end

for i=1:length(drug_concentrations)
    %if headers
    % Calculate APDs for this (last) stim.
    Upstroke(i) = data(start_idx+i-1, 3);
    Vmax(i) = data(start_idx+i-1, 4);
    APD50s(i) = data(start_idx+i-1, 5);
    APD90s(i) = data(start_idx+i-1, 6);
    AP_triangulation(i) = APD90s(i) - APD50s(i);
end

assert(length(APD90s)>=8)
assert(length(APD50s)>=8)

for i=1:7
    concs(i) = 10.0^(i-1);
end
figure(drug+15)
semilogx(concs,APD90s(2:8),'r.-',[1 1e7],[APD90s(1) APD90s(1)],'r--')
hold on
semilogx(concs,APD50s(2:8),'b.-',[1 1e7],[APD50s(1) APD50s(1)],'b--')
title(input_data.drug_names{drug})
legend('90','90 (no drug)','50','50 (no drug)')
xlim([1 1e6])
ylabel('Duration (ms)')
xlabel('Concentration')
