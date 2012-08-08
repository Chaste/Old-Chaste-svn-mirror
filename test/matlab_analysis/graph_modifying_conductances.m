function [ output_args ] = graph_modifying_conductances( input_args )
%GRAPH_MODIFYING_CONDUCTANCES Summary of this function goes here
%   Detailed explanation goes here
close all

data = load_drug_data();

data.num_stims = 1000; % Steady state now use 1000 stims.
do_S1S2_analysis = false;
do_dynamic_analysis = false;
plot_for_overdose_too = false;

base_results_folder = data.foldername;

num_models_to_analyse = [1 2 3 4 5];

% Load up pre-existing compiled_measures - we will overwrite
% whichever we are interested in below.
%load('classification/compiled_results.mat')

% IF YOU WANT TO SAVE COMPILED_RESULTS MAKE SURE IT IS NOT COMMENTED OUT
% AT THE BOTTOM OF THIS FILE...

for model_index = num_models_to_analyse % Loop over loading all model results to compile output for classification
    
    %% Simulation Parameters
    data.hertz = 1;
    % Get names and Column of data file corresponding to [Ca_i].
    % (Chaste ensures voltage is always the first variable in the output file)
    % Note that this is equal to C index + 2 as time column and indexing
    % from one instead of zero...

    if model_index==1  
        data.model_name = 'shannon_wang_puglisi_weber_bers_2004_model_updated'; % Full results
    elseif model_index==2
        data.model_name = 'tentusscher_panfilov_2006_epi'; % Full results
    elseif model_index==3
        data.model_name = 'hund_rudy_2004'; % Full results
    elseif model_index==4 
        data.model_name = 'mahajan_shiferaw_model_2008'; % Full results
    elseif model_index==5
        data.model_name = 'grandi2010'; % Full results
    else
        error('No such model')
    end
       
    %% Plotting Parameters
    plot_APD50_too = false;
    plot_HERG_only = false;
    reload_chaste_steady_state_data = false;
    reload_chaste_S1S2_data = false;
    reload_chaste_dynamic_data = false;
    
    if plot_HERG_only
        data.model_name = [data.model_name 'HergOnly'];
    end
    data.foldername = [base_results_folder data.model_name '/']; % Add model name to path
    
 % We only do the calculations for drugs with all three data points at the moment.        
    data.drug_indices_to_use = find((abs(data.drug_IC50_HERG + 1) > 1e-6)...
    & (abs(data.drug_IC50_Na + 1) > 1e-6) ...
    & (abs(data.drug_IC50_Ca + 1) > 1e-6));

    
%     if plot_HERG_only
%         % Figure out which drugs have HergIC50 data
%         % (not -1 in that column)
%         data.drug_indices_to_use = find(abs(data.drug_IC50_HERG + 1) > 1e-6)';
%     end
        
    %% Load Simulation Data
    if reload_chaste_steady_state_data
        % Plot all of the Action potentials for each concentration etc...
        for i = 1:length(data.drug_indices_to_use)
            drug = data.drug_indices_to_use(i);
            data.drug_names(drug)
            % Steady State Full Conductance Solution...
            [AP_triangulation{drug} APD90s{drug}, APD50s{drug} Vmax{drug} Upstroke{drug}] = LoadSteadyStateData(drug, data, 'voltage_results.dat');
            [Ca_triangulation{drug}, CaTransientD90s{drug}, CaTransientD50s{drug}, Ca_Peak_Info{drug}] = LoadSteadyStateData(drug, data, 'calcium_results.dat');
        end
        save(['results/' data.model_name '_' num2str(data.hertz) '_APD90s.mat'],'APD90s','-mat');
        save(['results/' data.model_name '_' num2str(data.hertz) '_APD50s.mat'],'APD50s','-mat');
        save(['results/' data.model_name '_' num2str(data.hertz) '_AP_triangulation.mat'],'AP_triangulation','-mat');
        save(['results/' data.model_name '_' num2str(data.hertz) '_Vmax.mat'],'Vmax','-mat');
        save(['results/' data.model_name '_' num2str(data.hertz) '_Upstroke.mat'],'Upstroke','-mat');
        save(['results/' data.model_name '_' num2str(data.hertz) '_CaTransientD90s.mat'],'CaTransientD90s','-mat');
        save(['results/' data.model_name '_' num2str(data.hertz) '_CaTransientD50s.mat'],'CaTransientD50s','-mat');
        save(['results/' data.model_name '_' num2str(data.hertz) '_Ca_triangulation.mat'],'Ca_triangulation','-mat');
        save(['results/' data.model_name '_' num2str(data.hertz) '_CaPeaks.mat'],'Ca_Peak_Info','-mat');
    else
        load(['results/' data.model_name '_' num2str(data.hertz) '_APD90s.mat']);
        load(['results/' data.model_name '_' num2str(data.hertz) '_APD50s.mat']);
        load(['results/' data.model_name '_' num2str(data.hertz) '_AP_triangulation.mat']);
        load(['results/' data.model_name '_' num2str(data.hertz) '_Vmax.mat']);
        load(['results/' data.model_name '_' num2str(data.hertz) '_Upstroke.mat']);
        load(['results/' data.model_name '_' num2str(data.hertz) '_CaTransientD90s.mat']);
        load(['results/' data.model_name '_' num2str(data.hertz) '_CaTransientD50s.mat']);
        load(['results/' data.model_name '_' num2str(data.hertz) '_Ca_triangulation.mat']);
        load(['results/' data.model_name '_' num2str(data.hertz) '_CaPeaks.mat']);
    end
    
    plot_steady_summary_graphs(data,[], APD90s, APD50s,'Voltage');
    
    %%
    if do_S1S2_analysis
        if reload_chaste_S1S2_data
            for i = 1:length(data.drug_indices_to_use)
                drug = data.drug_indices_to_use(i);
                data.drug_names(drug)
                % Load Chaste S1/S2 protocol results and save max restitution
                % slopes to file for this drug and concentration.
                [MaxSlopes{drug} S1S2APD90s{drug} S1S2PacingCycleLengths{drug}] = LoadS1S2data(drug, data);
            end
            save(['results/' data.model_name '_' num2str(data.hertz) '_MaxS1S2RestitutionSlope.mat'],'MaxSlopes','-mat');
            save(['results/' data.model_name '_' num2str(data.hertz) '_S1S2_APD90s.mat'],'S1S2APD90s','-mat');
            save(['results/' data.model_name '_' num2str(data.hertz) '_S1S2_PacingCycleLengths.mat'],'S1S2PacingCycleLengths','-mat');
        else
            load(['results/' data.model_name '_' num2str(data.hertz) '_MaxS1S2RestitutionSlope.mat']);
            load(['results/' data.model_name '_' num2str(data.hertz) '_S1S2_APD90s.mat']);
            load(['results/' data.model_name '_' num2str(data.hertz) '_S1S2_PacingCycleLengths.mat']);
        end
        %plot_s1s2_summary_graphs(data, S1S2APD90s, S1S2PacingCycleLengths, MaxSlopes);
    end
    
    %%
    if do_dynamic_analysis
        if reload_chaste_dynamic_data
            for i = 1:length(data.drug_indices_to_use)
                drug = data.drug_indices_to_use(i);
                data.drug_names(drug)
                % Load Chaste S1/S2 protocol results and save max restitution
                % slopes to file for this drug and concentration
                [DynamicMaxSlopes{drug} DynamicAlternansStartFrequency{drug} ...
                    DynamicEadOnsetFrequency{drug} DynamicAreas{drug} DynamicAPD90s{drug} DynamicPacingCycleLengths{drug}] = LoadDynamicRestitutionData(drug, data);
            end
            save(['results/' data.model_name '_' num2str(data.hertz) '_DynamicRestitutionMaxSlope.mat'],'DynamicMaxSlopes','-mat');
            save(['results/' data.model_name '_' num2str(data.hertz) '_DynamicAlternansStartFrequency.mat'],'DynamicAlternansStartFrequency','-mat');
            save(['results/' data.model_name '_' num2str(data.hertz) '_DynamicEadOnsetFrequency.mat'],'DynamicEadOnsetFrequency','-mat');
            save(['results/' data.model_name '_' num2str(data.hertz) '_DynamicAreas.mat'],'DynamicAreas','-mat');
            save(['results/' data.model_name '_' num2str(data.hertz) '_DynamicRestitutionAPD90s.mat'],'DynamicAPD90s','-mat');
            save(['results/' data.model_name '_' num2str(data.hertz) '_DynamicRestitutionPacingCycleLengths.mat'],'DynamicPacingCycleLengths','-mat');
        else
            load(['results/' data.model_name '_' num2str(data.hertz) '_DynamicRestitutionMaxSlope.mat']);
            load(['results/' data.model_name '_' num2str(data.hertz) '_DynamicAlternansStartFrequency.mat']);
            load(['results/' data.model_name '_' num2str(data.hertz) '_DynamicEadOnsetFrequency.mat']);
            load(['results/' data.model_name '_' num2str(data.hertz) '_DynamicAreas.mat']);
            load(['results/' data.model_name '_' num2str(data.hertz) '_DynamicRestitutionAPD90s.mat']);
            load(['results/' data.model_name '_' num2str(data.hertz) '_DynamicRestitutionPacingCycleLengths.mat']);
        end
        %plot_dynamic_summary_graphs(data, DynamicAPD90s, DynamicPacingCycleLengths);
    end
    
    %% Load Clinical Drug Data
    % Which drugs do we have clinical dose data for?
    drug_indices_with_clinical_doses = data.drug_indices_to_use(data.drug_low_ETPC(data.drug_indices_to_use)>0);
    plot_these_drugs = drug_indices_with_clinical_doses;
    add_to_figure_number = (model_index*11 + data.hertz*17);
    
    %% Plot the Raw IC50 values
    figure(120+add_to_figure_number)
    subplot(2,2,1)
    semilogy(data.drug_redfern(plot_these_drugs), data.drug_IC50_HERG(plot_these_drugs),'ko')
    ylabel('hERG IC_{50} (nM)','FontSize',16);
    subplot(2,2,2)
    semilogy(data.drug_redfern(plot_these_drugs), data.drug_IC50_Na(plot_these_drugs),'ko')
    ylabel('Na IC_{50} (nM)','FontSize',16);
    subplot(2,2,3)
    semilogy(data.drug_redfern(plot_these_drugs), data.drug_IC50_Ca(plot_these_drugs),'ko')
    ylabel('CaL IC_{50} (nM)','FontSize',16);
    subplot(2,2,4)
    semilogy(data.drug_redfern(plot_these_drugs), data.drug_high_ETPC(plot_these_drugs),'ko')
    ylabel('max EFTPC (nM)','FontSize',16);
    for i=1:4
        subplot(2,2,i)
        xlabel('Risk Category','FontSize',16);
        xlim([0.5 5.5])
        ylim([1e-1 10^6])
        set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
        set(gca,'YTick',[1e0 1e2 1e4 1e6],'FontSize',14);
    end
    
    % Plot the 3D representation of all three ic50s...
    figure(679)
    for i=1:length(data.drug_IC50_HERG)
        plot3(log(data.drug_IC50_HERG(i)),log(data.drug_IC50_Na(i)),...
            log(data.drug_IC50_Ca(i)), [data.drug_redfern_colour(i) 'o']);
        hold on
    end
    xlabel('log(Herg IC50)')
    ylabel('log(Na IC50)')
    zlabel('log(Ca IC50)')
    
%     % Plot the Lawrence markers...
%     figure(678)
%     for i=1:length(data.drug_IC50_HERG)
%         plot(data.drug_redfern, data.drug_triad_1ETPC,'bo',...
%             data.drug_redfern, data.drug_triad_10ETPC,'rx',...
%             data.drug_redfern, data.drug_triad_30ETPC,'gs',...
%             data.drug_redfern, data.drug_triad_100ETPC,'kd')
%     end
%     legend('EFTPC','10*EFTPC','30*EFTPC','100*EFTPC','Location','EastOutside')
%     xlabel('Risk category','FontSize',16)
%     ylabel('TRIad Score','FontSize',16)
%     xlim([0.5 5.5])
%     ylim([0 max(data.drug_triad_100ETPC)+5])
%     set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
%     title('SCREENIT score correlation with TdP risk')
    
    %% Plot the GSK Biomarkers
    figure(121+add_to_figure_number)
    subplot(1,2,1)
    semilogy(data.drug_redfern(plot_these_drugs), data.drug_IC50_Na_over_HERG(plot_these_drugs),'o')
    hold on
    ylabel('Na IC_{50}/HERG IC_{50}','FontSize',16);
    subplot(1,2,2)
    semilogy(data.drug_redfern(plot_these_drugs), data.drug_IC50_Ca_over_HERG(plot_these_drugs),'o')
    hold on
    ylabel('Ca IC_{50}/HERG IC_{50}','FontSize',16);
    for i=1:2
        subplot(1,2,i)
        xlabel('Risk category','FontSize',16);
        xlim([0.5 5.5])
        set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
        title('GSK Proposed Biomarker','FontSize',18)
        plot([0.5 5.5], [10 10],'k--')
    end
    
    %% Plot the APD values
    num_this_redfern = zeros(5,1);
    if plot_for_overdose_too
        num_doses = 4;
        dose_indices = 9:12;
    else
        num_doses = 3;
        dose_indices = 9:11;
    end
    
    %% Here we start to compile all of the possible measures into a single
    %%% structure for use in classification...
    %%% These measures are the same for every model
    compiled_results.measures(:,1) = data.drug_IC50_Na;
    compiled_results.measure_names{1} = 'Na IC50';
    compiled_results.measures(:,2) = data.drug_IC50_Ca;
    compiled_results.measure_names{2} = 'CaL IC50';
    compiled_results.measures(:,3) = data.drug_IC50_HERG;
    compiled_results.measure_names{3} = 'hERG IC50';
    compiled_results.measures(:,4) = data.drug_low_ETPC;
    compiled_results.measure_names{4} = 'EFTPC Low';
    compiled_results.measures(:,5) = data.drug_high_ETPC;
    compiled_results.measure_names{5} = 'EFTPC High';
    % The logs of the concentrations will probably be more interesting...
    compiled_results.measures(:,6) = log10(data.drug_IC50_Na);
    compiled_results.measure_names{6} = 'log10(Na IC50)';
    compiled_results.measures(:,7) = log10(data.drug_IC50_Ca);
    compiled_results.measure_names{7} = 'log10(CaL IC50)';
    compiled_results.measures(:,8) = log10(data.drug_IC50_HERG);
    compiled_results.measure_names{8} = 'log10(hERG IC50)';
    % Yi's measures:
    compiled_results.measures(:,9) = log10(data.drug_IC50_Na./data.drug_IC50_HERG);
    compiled_results.measure_names{9} = 'log10(Na IC50/hERG IC50)';
    compiled_results.measures(:,10) = log10(data.drug_IC50_Ca./data.drug_IC50_HERG);
    compiled_results.measure_names{10} = 'log10(Ca IC50/hERG IC50)';
    % Redfern Measure:
    compiled_results.measures(:,11) = log10(data.drug_IC50_HERG./data.drug_high_ETPC);
    compiled_results.measure_names{11} = 'log10(hERG IC50/high_ETPC)';
    
    %     compiled_results.measures(:,12) = data.drug_triad_1ETPC;
    %     compiled_results.measure_names{12} = 'Lawrence Triad at 1 EFTPC';
    %     compiled_results.measures(:,13) = data.drug_triad_10ETPC;
    %     compiled_results.measure_names{13} = 'Lawrence Triad at 10 EFTPC';
    %     compiled_results.measures(:,14) = data.drug_triad_30ETPC;
    %     compiled_results.measure_names{14} = 'Lawrence Triad at 30 EFTPC';
    %     compiled_results.measures(:,15) = data.drug_triad_100ETPC;
    %     compiled_results.measure_names{15} = 'Lawrence Triad at 100 EFTPC';
    
    %% These measures are model specific so add num_models + model index...
    num_measures = 75;
    start_index = 12 + (model_index-1)*num_measures;
    if plot_HERG_only
        start_index = start_index + length(num_models_to_analyse)*num_measures;
    end
    compiled_results.measure_names{start_index} = [data.model_name ' APD90 OverDose'];
    compiled_results.measure_names{start_index+1} = [data.model_name ' APD90 HighDose'];
    compiled_results.measure_names{start_index+2} = [data.model_name ' APD90 MedDose'];
    compiled_results.measure_names{start_index+3} = [data.model_name ' APD90 LowDose'];
    compiled_results.measure_names{start_index+4} = [data.model_name ' APD90 LargestEffectDose'];
    compiled_results.measure_names{start_index+5} = [data.model_name ' APD50 OverDose'];
    compiled_results.measure_names{start_index+6} = [data.model_name ' APD50 HighDose'];
    compiled_results.measure_names{start_index+7} = [data.model_name ' APD50 MedDose'];
    compiled_results.measure_names{start_index+8} = [data.model_name ' APD50 LowDose'];
    compiled_results.measure_names{start_index+9} = [data.model_name ' APD50 LargestEffectDose'];
    compiled_results.measure_names{start_index+10} = [data.model_name ' APD Triangulation OverDose'];
    compiled_results.measure_names{start_index+11} = [data.model_name ' APD Triangulation HighDose'];
    compiled_results.measure_names{start_index+12} = [data.model_name ' APD Triangulation MedDose'];
    compiled_results.measure_names{start_index+13} = [data.model_name ' APD Triangulation LowDose'];
    compiled_results.measure_names{start_index+14} = [data.model_name ' APD Triangulation LargestEffectDose'];
    compiled_results.measure_names{start_index+15} = [data.model_name ' Peak Vm OverDose'];
    compiled_results.measure_names{start_index+16} = [data.model_name ' Peak Vm HighDose'];
    compiled_results.measure_names{start_index+17} = [data.model_name ' Peak Vm MedDose'];
    compiled_results.measure_names{start_index+18} = [data.model_name ' Peak Vm LowDose'];
    compiled_results.measure_names{start_index+19} = [data.model_name ' Peak Vm LargestEffectDose'];
    compiled_results.measure_names{start_index+20} = [data.model_name ' Upstoke Velocity OverDose'];
    compiled_results.measure_names{start_index+21} = [data.model_name ' Upstoke Velocity HighDose'];
    compiled_results.measure_names{start_index+22} = [data.model_name ' Upstoke Velocity MedDose'];
    compiled_results.measure_names{start_index+23} = [data.model_name ' Upstoke Velocity LowDose'];
    compiled_results.measure_names{start_index+24} = [data.model_name ' Upstoke Velocity LargestEffectDose'];
    compiled_results.measure_names{start_index+25} = [data.model_name ' CaD90 OverDose'];
    compiled_results.measure_names{start_index+26} = [data.model_name ' CaD90 HighDose'];
    compiled_results.measure_names{start_index+27} = [data.model_name ' CaD90 MedDose'];
    compiled_results.measure_names{start_index+28} = [data.model_name ' CaD90 LowDose'];
    compiled_results.measure_names{start_index+29} = [data.model_name ' CaD90 LargestEffectDose'];
    compiled_results.measure_names{start_index+30} = [data.model_name ' CaD50 OverDose'];
    compiled_results.measure_names{start_index+31} = [data.model_name ' CaD50 HighDose'];
    compiled_results.measure_names{start_index+32} = [data.model_name ' CaD50 MedDose'];
    compiled_results.measure_names{start_index+33} = [data.model_name ' CaD50 LowDose'];
    compiled_results.measure_names{start_index+34} = [data.model_name ' CaD50 LargestEffectDose'];
    compiled_results.measure_names{start_index+35} = [data.model_name ' CaD Triangulation OverDose'];
    compiled_results.measure_names{start_index+36} = [data.model_name ' CaD Triangulation HighDose'];
    compiled_results.measure_names{start_index+37} = [data.model_name ' CaD Triangulation MedDose'];
    compiled_results.measure_names{start_index+38} = [data.model_name ' CaD Triangulation LowDose'];
    compiled_results.measure_names{start_index+39} = [data.model_name ' CaD Triangulation LargestEffectDose'];
    compiled_results.measure_names{start_index+40} = [data.model_name ' Peak Ca OverDose'];
    compiled_results.measure_names{start_index+41} = [data.model_name ' Peak Ca HighDose'];
    compiled_results.measure_names{start_index+42} = [data.model_name ' Peak Ca MedDose'];
    compiled_results.measure_names{start_index+43} = [data.model_name ' Peak Ca LowDose'];
    compiled_results.measure_names{start_index+44} = [data.model_name ' Peak Ca LargestEffectDose'];    
    compiled_results.measure_names{start_index+45} = [data.model_name ' S1S2 Max Slope OverDose'];
    compiled_results.measure_names{start_index+46} = [data.model_name ' S1S2 Max Slope HighDose'];
    compiled_results.measure_names{start_index+47} = [data.model_name ' S1S2 Max Slope MedDose'];
    compiled_results.measure_names{start_index+48} = [data.model_name ' S1S2 Max Slope LowDose'];
    compiled_results.measure_names{start_index+49} = [data.model_name ' S1S2 Max Slope LargestEffectDose'];
    compiled_results.measure_names{start_index+50} = [data.model_name ' Dynamic Max Slope OverDose'];
    compiled_results.measure_names{start_index+51} = [data.model_name ' Dynamic Max Slope HighDose'];
    compiled_results.measure_names{start_index+52} = [data.model_name ' Dynamic Max Slope MedDose'];
    compiled_results.measure_names{start_index+53} = [data.model_name ' Dynamic Max Slope LowDose'];
    compiled_results.measure_names{start_index+54} = [data.model_name ' Dynamic Max Slope LargestEffectDose'];
    compiled_results.measure_names{start_index+55} = [data.model_name ' Dynamic Alternans Start Freq OverDose'];
    compiled_results.measure_names{start_index+56} = [data.model_name ' Dynamic Alternans Start Freq HighDose'];
    compiled_results.measure_names{start_index+57} = [data.model_name ' Dynamic Alternans Start Freq MedDose'];
    compiled_results.measure_names{start_index+58} = [data.model_name ' Dynamic Alternans Start Freq LowDose'];
    compiled_results.measure_names{start_index+59} = [data.model_name ' Dynamic Alternans Start Freq LargestEffectDose'];
    compiled_results.measure_names{start_index+60} = [data.model_name ' Dynamic EAD Start Freq OverDose'];
    compiled_results.measure_names{start_index+61} = [data.model_name ' Dynamic EAD Start Freq HighDose'];
    compiled_results.measure_names{start_index+62} = [data.model_name ' Dynamic EAD Start Freq MedDose'];
    compiled_results.measure_names{start_index+63} = [data.model_name ' Dynamic EAD Start Freq LowDose'];
    compiled_results.measure_names{start_index+64} = [data.model_name ' Dynamic EAD Start Freq LargestEffectDose'];
    compiled_results.measure_names{start_index+65} = [data.model_name ' Dynamic Instability Onset Freq OverDose'];
    compiled_results.measure_names{start_index+66} = [data.model_name ' Dynamic Instability Onset Freq HighDose'];
    compiled_results.measure_names{start_index+67} = [data.model_name ' Dynamic Instability Onset Freq MedDose'];
    compiled_results.measure_names{start_index+68} = [data.model_name ' Dynamic Instability Onset Freq LowDose'];
    compiled_results.measure_names{start_index+69} = [data.model_name ' Dynamic Instability Onset Freq LargestEffectDose'];
    compiled_results.measure_names{start_index+70} = [data.model_name ' Dynamic Area Between Curves OverDose'];
    compiled_results.measure_names{start_index+71} = [data.model_name ' Dynamic Area Between Curves HighDose'];
    compiled_results.measure_names{start_index+72} = [data.model_name ' Dynamic Area Between Curves MedDose'];
    compiled_results.measure_names{start_index+73} = [data.model_name ' Dynamic Area Between Curves LowDose'];
    compiled_results.measure_names{start_index+74} = [data.model_name ' Dynamic Area Between Curves LargestEffectDose'];
    
    %%
    for j=1:length(plot_these_drugs)
        % Get the global drug index
        drug = plot_these_drugs(j);
        clinical_concs(1) = data.drug_low_ETPC(drug);
        clinical_concs(2) = (data.drug_low_ETPC(drug)+data.drug_high_ETPC(drug))/2.0;
        clinical_concs(3) = data.drug_high_ETPC(drug);
        if plot_for_overdose_too
            clinical_concs(4) = 10*data.drug_high_ETPC(drug);
        end
        % This is just to make the plot look pretty with no overlaps
        redfern = data.drug_redfern(drug);
        ref_x = redfern-0.15;
        num_this_redfern(redfern) = num_this_redfern(redfern)+1;
        redfern_x = ref_x + 0.05*num_this_redfern(redfern);
        
        % Do the Redfern 30-fold plot
        figure(122+add_to_figure_number)
        plot_redfern_style(redfern_x,log10(data.drug_IC50_HERG(drug)./clinical_concs))
        
        % Retrieve clinical APDs
        APD90 = APD90s{drug};
        APD50 = APD50s{drug};
        APTriangulation = AP_triangulation{drug};
        PeakVm = Vmax{drug};
        UpstrokeVelocity = Upstroke{drug};
        CaTransientD90 = CaTransientD90s{drug};
        CaTransientD50 = CaTransientD50s{drug};
        CaTriangulation = Ca_triangulation{drug};
        CaPeaks = Ca_Peak_Info{drug};
        assert(length(APD90)==12)
        assert(length(APD50)==12)
        assert(length(CaTransientD90)==12)
        assert(length(CaTransientD50)==12)
        
        figure(123+add_to_figure_number)
        subplot(1,2,1)
        plot_redfern_style(redfern_x, APD90(dose_indices))
        if plot_APD50_too
            plot_redfern_style(redfern_x,APD50(dose_indices))
        end
        subplot(1,2,2)
        plot_redfern_style(redfern_x, 100*APD90(dose_indices)./APD90(1))
        if plot_APD50_too
            plot_redfern_style(redfern_x, 100*APD50(dose_indices)./APD50(1))
        end
        
        figure(139+add_to_figure_number)
        subplot(1,2,1)
        plot_redfern_style(redfern_x, PeakVm(dose_indices))
        subplot(1,2,2)
        plot_redfern_style(redfern_x, UpstrokeVelocity(dose_indices))
        
        compiled_results.measures(drug,start_index) = APD90(12); % overdose
        compiled_results.measures(drug,start_index + 1) = APD90(11);
        compiled_results.measures(drug,start_index + 2) = APD90(10);
        compiled_results.measures(drug,start_index + 3) = APD90(9);
        compiled_results.measures(drug,start_index + 4) = max([APD90(9)-APD90(1) APD90(10)-APD90(1) APD90(11)-APD90(1)]);
        compiled_results.measures(drug,start_index + 5) = APD50(12);
        compiled_results.measures(drug,start_index + 6) = APD50(11);
        compiled_results.measures(drug,start_index + 7) = APD50(10);
        compiled_results.measures(drug,start_index + 8) = APD50(9);
        compiled_results.measures(drug,start_index + 9) = max([APD50(9)-APD50(1) APD50(10)-APD50(1) APD50(11)-APD50(1)]);
        compiled_results.measures(drug,start_index + 10) = APTriangulation(12);
        compiled_results.measures(drug,start_index + 11) = APTriangulation(11);
        compiled_results.measures(drug,start_index + 12) = APTriangulation(10);
        compiled_results.measures(drug,start_index + 13) = APTriangulation(9);
        [~,temp] = max([abs(APTriangulation(9)-APTriangulation(1)) abs(APTriangulation(10)-APTriangulation(1)) abs(APTriangulation(11)-APTriangulation(1))]);
        compiled_results.measures(drug,start_index + 14) = APTriangulation(8+temp);
        compiled_results.measures(drug,start_index + 15) = PeakVm(12);
        compiled_results.measures(drug,start_index + 16) = PeakVm(11);
        compiled_results.measures(drug,start_index + 17) = PeakVm(10);
        compiled_results.measures(drug,start_index + 18) = PeakVm(9);
        [~,temp] = max([abs(PeakVm(9)-PeakVm(1)) abs(PeakVm(10)-PeakVm(1)) abs(PeakVm(11)-PeakVm(1))]);
        compiled_results.measures(drug,start_index + 19) = PeakVm(8+temp);
        compiled_results.measures(drug,start_index + 20) = UpstrokeVelocity(12);
        compiled_results.measures(drug,start_index + 21) = UpstrokeVelocity(11);
        compiled_results.measures(drug,start_index + 22) = UpstrokeVelocity(10);
        compiled_results.measures(drug,start_index + 23) = UpstrokeVelocity(9);
        [~,temp] = max([abs(UpstrokeVelocity(9)-UpstrokeVelocity(1)) abs(UpstrokeVelocity(10)-UpstrokeVelocity(1)) abs(UpstrokeVelocity(11)-UpstrokeVelocity(1)) ]);
        compiled_results.measures(drug,start_index + 24) = UpstrokeVelocity(8+temp);
        compiled_results.measures(drug,start_index + 25) = CaTransientD90(12); % overdose
        compiled_results.measures(drug,start_index + 26) = CaTransientD90(11);
        compiled_results.measures(drug,start_index + 27) = CaTransientD90(10);
        compiled_results.measures(drug,start_index + 28) = CaTransientD90(9);
        compiled_results.measures(drug,start_index + 29) = max([CaTransientD90(9)-CaTransientD90(1) CaTransientD90(10)-CaTransientD90(1) CaTransientD90(11)-CaTransientD90(1)]);
        compiled_results.measures(drug,start_index + 30) = CaTransientD50(12);
        compiled_results.measures(drug,start_index + 31) = CaTransientD50(11);
        compiled_results.measures(drug,start_index + 32) = CaTransientD50(10);
        compiled_results.measures(drug,start_index + 33) = CaTransientD50(9);
        compiled_results.measures(drug,start_index + 34) = max([CaTransientD50(9)-CaTransientD50(1) CaTransientD50(10)-CaTransientD50(1) CaTransientD50(11)-CaTransientD50(1)]);
        compiled_results.measures(drug,start_index + 35) = CaTriangulation(12);
        compiled_results.measures(drug,start_index + 36) = CaTriangulation(11);
        compiled_results.measures(drug,start_index + 37) = CaTriangulation(10);
        compiled_results.measures(drug,start_index + 38) = CaTriangulation(9);
        [~,temp] = max([abs(CaTriangulation(9)-CaTriangulation(1)) abs(CaTriangulation(10)-CaTriangulation(1)) abs(CaTriangulation(11)-CaTriangulation(1))]);
        compiled_results.measures(drug,start_index + 39) = CaTriangulation(8+temp);
        compiled_results.measures(drug,start_index + 40) = CaPeaks(12);
        compiled_results.measures(drug,start_index + 41) = CaPeaks(11);
        compiled_results.measures(drug,start_index + 42) = CaPeaks(10);
        compiled_results.measures(drug,start_index + 43) = CaPeaks(9);
        [~,temp] = max([abs(CaPeaks(9)-CaPeaks(1)) abs(CaPeaks(10)-CaPeaks(1)) abs(CaPeaks(11)-CaPeaks(1))]);
        compiled_results.measures(drug,start_index + 44) = CaPeaks(8+temp);
        
        figure(678)
        plot(ones(num_doses,1).*log10(data.drug_IC50_HERG(drug)),100*APD90(dose_indices)./APD90(1), [data.drug_redfern_colour(drug) '.-']);
        xlabel('log10(Herg IC50)')
        ylabel('% change in APD')
        hold on
        
        % Calcium Transient Duration plots
        figure(136+add_to_figure_number)
        subplot(1,2,1)
        plot_redfern_style(redfern_x,CaTransientD90(dose_indices))
        if plot_APD50_too
            plot_redfern_style(redfern_x,CaTransientD50(dose_indices))
        end
        subplot(1,2,2)
        plot_redfern_style(redfern_x,100*CaTransientD90(dose_indices)./CaTransientD90(1)); % Just plot clinical APDs
        if plot_APD50_too
            plot_redfern_style(redfern_x,100*CaTransientD50(dose_indices)./CaTransientD50(1)); % Just plot clinical APDs
        end

        % Calcium Transient Peak plots
        figure(137+add_to_figure_number)
        subplot(1,2,1)
        plot_redfern_style(redfern_x,CaPeaks(dose_indices)) % Just plot clinical APDs
        subplot(1,2,2)
        plot_redfern_style(redfern_x,100*CaPeaks(dose_indices)./CaPeaks(1)); % Just plot clinical APDs
        
        
        % Calcium Triangulation plots
        figure(138+add_to_figure_number)
        subplot(1,2,1)
        plot_redfern_style(redfern_x,CaTriangulation(dose_indices))
        subplot(1,2,2)
        plot_redfern_style(redfern_x,100*CaTriangulation(dose_indices)./CaTriangulation(1))% Just plot clinical APDs
        
        %     % Special plot for results collection.
        %     %hold on
        %     if find(drug_indices_with_multichannel_data==drug)
        %         colour = 'r.-';
        %     else
        %         colour = 'k.-';
        %     end
        %     figure(a)
        %     plot(ones(num_doses,1).*redfern_x, 100*APD90(dose_indices)./APD90(1),colour);
        %     hold on
        %     xlabel('Risk category','FontSize',16);
        %     xlim([0.5 5.5])
        %     set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
        %     title_text = [data.model_name ' APD90 Change'];
        %     if plot_HERG_only
        %         title_text = [title_text ' (Herg-only block)'];
        %     else
        %         title_text = [title_text ' (Multichannel block)'];
        %     end
        %     title(title_text,'FontSize',18)
        %     ylabel('Steady State APD90 (% of drug free APD)','FontSize',16);
        %     plot([0.5 5.5], [100 100],'k--')
        
        %% AP triangulation
        figure(126+add_to_figure_number)
        subplot(1,2,1)
        plot_redfern_style(redfern_x,APTriangulation(dose_indices))
        subplot(1,2,2)
        plot_redfern_style(redfern_x,100.0*APTriangulation(dose_indices)./APTriangulation(1))
        
        %     figure(130+add_to_figure_number)
        %     plot3(ones(8,1).*x,log(data.drug_concentrations),100*APD90(1:8)./APD90(1),'b.-')
        %     hold on
        %     plot3(ones(3,1).*x,log(clinical_concs),100*APD50(9:11)./APD50(1),'r.-')
        %
        %     figure(131+add_to_figure_number)
        %     red_val = (redfern-0.05 + 0.05*num_this_redfern);
        %     surfc([red_val-0.025 red_val+0.025],log(clinical_concs),[100*APD90(9:11)./APD90(1); 100*APD90(9:11)./APD90(1)]')
        %     hold on
        
        if do_S1S2_analysis
            MaxSlope = MaxSlopes{drug};
            assert(length(MaxSlope)==12)
            figure(124+add_to_figure_number)
            % For all drug concentrations...
            %plot(ones(8,1).*redfern_x, MaxSlope(1:8),'b.-');
            %hold on
            plot_redfern_style(redfern_x,MaxSlope(dose_indices))
            compiled_results.measures(drug,start_index + 45) = MaxSlope(12);
            compiled_results.measures(drug,start_index + 46) = MaxSlope(11);
            compiled_results.measures(drug,start_index + 47) = MaxSlope(10);
            compiled_results.measures(drug,start_index + 48) = MaxSlope(9);
            [~,temp] = max([abs(MaxSlope(9)-MaxSlope(1)) abs(MaxSlope(10)-MaxSlope(1)) abs(MaxSlope(11)-MaxSlope(1)) ]);
            compiled_results.measures(drug,start_index + 49) = MaxSlope(8+temp);
        end
        
        
        if do_dynamic_analysis
            MaxSlope = DynamicMaxSlopes{drug};
            assert(length(MaxSlope)==12)
            figure(125+add_to_figure_number)
            % For all drug concentrations...
            %plot(ones(8,1).*redfern_x, MaxSlope(1:8),'b.-');
            %hold on
            plot_redfern_style(redfern_x,MaxSlope(dose_indices))
            compiled_results.measures(drug,start_index + 50) = MaxSlope(12);
            compiled_results.measures(drug,start_index + 51) = MaxSlope(11);
            compiled_results.measures(drug,start_index + 52) = MaxSlope(10);
            compiled_results.measures(drug,start_index + 53) = MaxSlope(9);
            [~,temp] = max([abs(MaxSlope(9)-MaxSlope(1)) abs(MaxSlope(10)-MaxSlope(1)) abs(MaxSlope(11)-MaxSlope(1)) ]);
            compiled_results.measures(drug,start_index + 54) = MaxSlope(8+temp);
            
            alternans_start = DynamicAlternansStartFrequency{drug};
            
            figure(127+add_to_figure_number)
            alternans_starts = [alternans_start{9} alternans_start{10} alternans_start{11}];
            plot_redfern_style(redfern_x, alternans_starts)
            alternans_starts = [alternans_starts alternans_start{12}];
            if plot_for_overdose_too % Replot with extended data...
                plot_redfern_style(redfern_x, alternans_starts)
            end
            compiled_results.measures(drug,start_index + 55) = alternans_starts(4);
            compiled_results.measures(drug,start_index + 56) = alternans_starts(3);
            compiled_results.measures(drug,start_index + 57) = alternans_starts(2);
            compiled_results.measures(drug,start_index + 58) = alternans_starts(1);
            [~,temp] = max([alternans_starts(1) alternans_starts(2) alternans_starts(3)]);
            compiled_results.measures(drug,start_index + 59) = alternans_starts(temp);
            
            
            figure(128+add_to_figure_number)
            ead_start = DynamicEadOnsetFrequency{drug};
            ead_starts = [ead_start{9} ead_start{10} ead_start{11} ];
            plot_redfern_style(redfern_x, ead_starts)
            ead_starts = [ead_starts ead_start{12}];
            if plot_for_overdose_too
                plot_redfern_style(redfern_x, ead_starts)
            end                      
            compiled_results.measures(drug,start_index + 60) = ead_starts(4);
            compiled_results.measures(drug,start_index + 61) = ead_starts(3);
            compiled_results.measures(drug,start_index + 62) = ead_starts(2);
            compiled_results.measures(drug,start_index + 63) = ead_starts(1);
            [~,temp] = max([ead_starts(1) ead_starts(2) ead_starts(3)]);           
            compiled_results.measures(drug,start_index + 64) = ead_starts(temp);
            
            for i=1:length(alternans_starts)
                onset_of_instability(i) = max([alternans_starts(i) ead_starts(i)]);
            end
            
            figure(129+add_to_figure_number)
            plot_redfern_style(redfern_x,onset_of_instability(1:num_doses))
            compiled_results.measures(drug,start_index + 65) = onset_of_instability(4);
            compiled_results.measures(drug,start_index + 66) = onset_of_instability(3);
            compiled_results.measures(drug,start_index + 67) = onset_of_instability(2);
            compiled_results.measures(drug,start_index + 68) = onset_of_instability(1);
            [~,temp] = max([onset_of_instability(1) onset_of_instability(2) onset_of_instability(3)]);       
            compiled_results.measures(drug,start_index + 69) = onset_of_instability(temp);
            
            areas = DynamicAreas{drug};
            figure(132+add_to_figure_number)
            plot_redfern_style(redfern_x,areas(dose_indices))
            compiled_results.measures(drug,start_index + 70) = areas(12);
            compiled_results.measures(drug,start_index + 71) = areas(11);
            compiled_results.measures(drug,start_index + 72) = areas(10);
            compiled_results.measures(drug,start_index + 73) = areas(9);
            [~,temp] = max([areas(9)-areas(1) areas(10)-areas(1) areas(11)-areas(1) ]);
            compiled_results.measures(drug,start_index + 74) = areas(8+temp);
        end
    end
    
    figure(122+add_to_figure_number)
    plot([0.5 5.5],log10([30 30]),'k--')
    xlim([0.5 5.5])
    set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
    xlabel('Risk category','FontSize',16);
    ylabel('log10(IC50/[EFTPC])','FontSize',16);
    title('Redfern Safety Measure','FontSize',18)
    legend(data.drug_names(drug_indices_with_clinical_doses),'Location','EastOutside');
    
    % Print axes etc.
    for j=[123 136 138]
        figure(j+add_to_figure_number)
        if j==123
            textdesc = 'APD90';
        elseif j==136
            textdesc = 'Ca transient D90';
        elseif j==138
            textdesc = 'Ca Triangulation';
        end
        for i=1:2
            subplot(1,2,i)
            xlabel('Risk category','FontSize',16);
            xlim([0.5 5.5])
            set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
            if i==1
                title([data.model_name ' ' num2str(data.hertz) 'Hz simulation of ' textdesc],'FontSize',18)
                ylabel(['Steady State ' textdesc ' under clinical drug concentrations'],'FontSize',16);
            elseif i==2
                title_text = [textdesc ' Change'];
                if plot_HERG_only
                    title_text = [title_text ' (Herg-only block)'];
                else
                    title_text = [title_text ' (Multichannel block)'];
                end
                title(title_text,'FontSize',18)
                ylabel(['Steady State ' textdesc ' (% of drug free ' textdesc ')'],'FontSize',16);
                plot([0.5 5.5], [100 100],'k--')
                legend(data.drug_names(drug_indices_with_clinical_doses),'Location','EastOutside')
            end
        end
    end
    
    figure(137+add_to_figure_number)
    for i=1:2
        subplot(1,2,i)
        xlabel('Risk category','FontSize',16);
        xlim([0.5 5.5])
        set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
        if i==1
            title([data.model_name ' ' num2str(data.hertz) 'Hz simulation of peak [Ca]i'],'FontSize',18)
            ylabel('Peak value of intracellular Calcium','FontSize',16);
            plot([0.5 5.5],[CaPeaks(1) CaPeaks(1)],'k--')
        elseif i==2
            title_text = ['Time of Peak Calcium [Ca]i'];
            if plot_HERG_only
                title_text = [title_text ' (Herg-only block)'];
            else
                title_text = [title_text ' (Multichannel block)'];
            end
            title(title_text,'FontSize',18)
            ylabel('Time of peak Calcium (ms)','FontSize',16);
            legend(data.drug_names(drug_indices_with_clinical_doses),'Location','EastOutside')
        end
    end
    
    figure(126+add_to_figure_number)
    for i=1:2
        subplot(1,2,i)
        xlabel('Risk category','FontSize',16);
        xlim([0.5 5.5])
        set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
        if i==1
            title([data.model_name ' ' num2str(data.hertz) 'Hz, APD Triangulation'],'FontSize',18)
            ylabel('APD90-APD50','FontSize',16);
        else
            title_text = 'Relative APD Triangulation';
            if plot_HERG_only
                title_text = [title_text ' (Herg-only block)'];
            else
                title_text = [title_text ' (Multichannel block)'];
            end
            title(title_text,'FontSize',18)
            ylabel('APD90-APD50 (% of drug free APD90-APD50)','FontSize',16);
            plot([0.5 5.5], [100 100],'k--')
            legend(data.drug_names(drug_indices_with_clinical_doses),'Location','EastOutside')
        end
    end
    
    figure(139+add_to_figure_number)
    for i=1:2
        subplot(1,2,i)
        xlabel('Risk category','FontSize',16);
        xlim([0.5 5.5])
        set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
        if i==1
            title([data.model_name ' ' num2str(data.hertz) 'Hz, Peak Vm'],'FontSize',18)
            ylabel('Peak Vm (mV)','FontSize',16);
        else
            title([data.model_name ' ' num2str(data.hertz) 'Hz, Max Upstroke Velocity (mV/ms)'],'FontSize',18)
            ylabel('Max Upstroke Velocity (mV/ms)','FontSize',16);
            plot([0.5 5.5], [100 100],'k--')
            legend(data.drug_names(drug_indices_with_clinical_doses),'Location','EastOutside')
        end
    end
    
    % for i=0:1
    % figure(130+i+add_to_figure_number)
    % xlabel('Risk category','FontSize',16);
    % xlim([0.5 5.5])
    % set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
    % ylabel('log(Drug Conc)','FontSize',16)
    % zlabel('Steady State APD (% of drug free APD)','FontSize',16);
    % end
    
    if do_S1S2_analysis
        % Print axes etc.
        figure(124+add_to_figure_number)
        
        xlabel('Risk category','FontSize',16);
        xlim([0.5 5.5])
        set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
        title([data.model_name ' with S1 = ' num2str(data.hertz) 'Hz'],'FontSize',18)
        ylabel('Maximum Slope of S1-S2 Restitution Curve','FontSize',16);
        legend(data.drug_names(drug_indices_with_clinical_doses),'Location','EastOutside')
    end
    
    
    if do_dynamic_analysis
        % Print axes etc.
        figure(125+add_to_figure_number)
        xlabel('Risk category','FontSize',16);
        xlim([0.5 5.5]);
        set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
        title([data.model_name ' simulation of APD'],'FontSize',18)
        ylabel('Maximum Slope of Dynamic Restitution Curve','FontSize',16);
        legend(data.drug_names(drug_indices_with_clinical_doses),'Location','EastOutside')
        
        
        for i=1:3
            figure(126+i+add_to_figure_number)
            xlabel('Risk category','FontSize',16);
            xlim([0.5 5.5]);
            set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
            current_ylim = get(gca,'ylim');
            ylim([current_ylim(1)-10 current_ylim(2)+25]);
            title([data.model_name ' dynamic restitution simulation'],'FontSize',18)
            set(gca,'YGrid','on','YTick',sort(DynamicPacingCycleLengths{1}));
            if i==1
                ylabel('Alternans onset pacing period (ms)','FontSize',16);
            elseif i==2
                ylabel('EAD onset pacing period (ms)','FontSize',16);
            elseif i==3
                ylabel('Instability onset pacing period (ms)','FontSize',16);
            end
            legend(data.drug_names(drug_indices_with_clinical_doses),'Location','EastOutside');
        end
        
        figure(132+add_to_figure_number)
        xlabel('Risk category','FontSize',16);
        xlim([0.5 5.5]);
        set(gca,'XTick',[1 2 3 4 5],'FontSize',14);
        title([data.model_name ' simulation of APD'],'FontSize',18)
        ylabel('Area between dynamic restitution curves (control and drug)','FontSize',16);
        legend(data.drug_names(drug_indices_with_clinical_doses),'Location','EastOutside')
        plot([-1 6],[0 0],'k--')
        text = data.drug_names(drug_indices_with_clinical_doses)
    end
end

%compiled_results.measure_names'
%length(compiled_results.measure_names)
%save('classification/compiled_results.mat','compiled_results','-mat')


function plot_redfern_style(redfern_x, dataset)
num_doses = length(dataset);
plot(ones(num_doses,1).*redfern_x, dataset,'k-'); % Just plot clinical APDs
hold all

plot(redfern_x, dataset(1),'ko','MarkerSize',1.0,'MarkerFaceColor','k'); % Just plot clinical APDs
plot(redfern_x, dataset(2),'ko','MarkerSize',3.0,'MarkerFaceColor','k'); % Just plot clinical APDs
plot(redfern_x, dataset(3),'ko','MarkerSize',4.0,'MarkerFaceColor','k'); % Just plot clinical APDs


