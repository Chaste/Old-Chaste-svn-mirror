function [AP_triangulation controlAPD90_50 peak_info] = plot_steady_summary_graphs(data, All_AP_traces, APD90s, APD50s, VoltageOrCalcium)
concs(1) = 0.0;
for i=2:8
    concs(i) = 10.0^(i-2);
end
last_entry = data.drug_indices_to_use(end);
size_of_plots = ceil(sqrt(last_entry));

if (strcmp(VoltageOrCalcium,'Voltage')==1)
    title_text = 'APD';
    figure_add = 0;
elseif (strcmp(VoltageOrCalcium,'Calcium')==1)
    title_text = 'Ca dur';
    figure_add = 400;
else
    error('Unrecognised option')
end

for i = 1:length(data.drug_indices_to_use)
    drug = data.drug_indices_to_use(i);
    controlAPD90_50 = APD90s{drug}-APD50s{drug};
    controlAPD90_50 = controlAPD90_50(1); % First entry is drug free...
    concentration = data.drug_concentrations_text;
    
    extra_concs_low = data.drug_low_ETPC;
    extra_concs_high = data.drug_high_ETPC;
    if extra_concs_low(drug) > 0 % If data is missing it is given a '-1' in the file
        more_concs(1) = extra_concs_low(drug);
        more_concs(2) = (extra_concs_low(drug)+extra_concs_high(drug))/2.0;
        more_concs(3) = extra_concs_high(drug);
        more_concs(4) = 10*extra_concs_high(drug);
        concentration{length(data.drug_concentrations_text)+1} = num2str(more_concs(1));
        concentration{length(data.drug_concentrations_text)+2} = num2str(more_concs(2));
        concentration{length(data.drug_concentrations_text)+3} = num2str(more_concs(3));
        concentration{length(data.drug_concentrations_text)+4} = num2str(more_concs(4));
    end

    fignum = ceil(drug/4);
    num_cols_on_big_fig = 4;
       
    basic_conc_indices = 1:length(concs);
    clinical_conc_indices = length(concs)+1:length(concs)+4;

    figure(120 + figure_add)
    subplot(size_of_plots,size_of_plots,drug)
    apd90s = APD90s{drug};
    apd50s = APD50s{drug};
    semilogx(concs, apd90s(basic_conc_indices),'b.-')
    hold on
    semilogx(concs, apd50s(basic_conc_indices),'k.-')
    if extra_concs_low(drug) > 0
        semilogx(more_concs, apd90s(clinical_conc_indices),'r.-')
        semilogx(more_concs, apd50s(clinical_conc_indices),'r.-')
    end
    ylabel(title_text,'FontSize',14)
    xlabel('[Drug]','FontSize',14)
    xlim([1 1e6])
    title([data.drug_names{drug}],'FontSize',16)
    
    figure(121 + figure_add)
    subplot(size_of_plots,size_of_plots,drug)
    apd90s = APD90s{drug};
    apd50s = APD50s{drug};
    control90 = apd90s(1);
    control50 = apd50s(1);
    semilogx(concs, 100*apd90s(basic_conc_indices)./control90,'b.-')
    hold on
    semilogx(concs, 100*apd50s(basic_conc_indices)./control50,'k.-')
    if extra_concs_low(drug) > 0
        semilogx(more_concs, 100*apd90s(clinical_conc_indices)./control90,'r.-')
        semilogx(more_concs, 100*apd50s(clinical_conc_indices)./control50,'r.-')
    end
    ylabel([title_text ' (% control)'],'FontSize',14)
    xlabel('[Drug]','FontSize',14)
    xlim([1 1e6])
    title([data.drug_names{drug}],'FontSize',16)
    
    if length(apd90s)==8
        apd90s = [apd90s 0 0 0 0];
    end    
    if length(apd50s)==8
        apd50s = [apd50s 0 0 0 0];
    end
    assert(length(apd90s)==12)
    assert(length(apd50s)==12)
    apd50s_for_excel(:,drug) = apd50s;
    apd90s_for_excel(:,drug) = apd90s;

       
    figure(122 + figure_add)
    subplot(size_of_plots,size_of_plots,drug)
    semilogx(concs, apd90s(basic_conc_indices)-apd50s(basic_conc_indices),'b.-')
    if extra_concs_low(drug) > 0
        hold on
        semilogx(more_concs, apd90s(clinical_conc_indices)-apd50s(clinical_conc_indices),'r.-')
    end
    ylabel([title_text ' 90-50'],'FontSize',14)
    xlabel('[Drug]','FontSize',14)
    xlim([1 1e6])
    title([data.drug_names{drug}],'FontSize',16)

    AP_triangulation{drug}=[];
    figure(123 + figure_add)
    subplot(size_of_plots,size_of_plots,drug)
    semilogx(concs, (apd90s(basic_conc_indices)-apd50s(basic_conc_indices))./controlAPD90_50,'b.-')
    if extra_concs_low(drug) > 0
        hold on
        semilogx(more_concs, (apd90s(clinical_conc_indices)-apd50s(clinical_conc_indices))./controlAPD90_50,'r.-')
        AP_triangulation{drug} = apd90s(clinical_conc_indices)-apd50s(clinical_conc_indices);
    end
    xlabel('[Drug]','FontSize',14)
    xlim([1 1e6])
    ylabel(['(' title_text ' 90-50)/control'],'FontSize',14)
    title([data.drug_names{drug}])
end

%apd50s_for_excel

