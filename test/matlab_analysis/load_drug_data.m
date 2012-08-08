function data = load_drug_data( plot_HERG_only )
%LOAD_DRUG_DATA Loads the drug data from the file drug_data.dat
% that contains the IC50 values and Redfern categories for each drug.

% Import the file
fid = fopen('../drug_data.dat');
C = textscan(fid, '%s %f %f %f %f %f %f %s %f ');
fclose(fid);

data.foldername = [getenv('CHASTE_TEST_OUTPUT') '/VaryingConductances/'];

data.drug_names = C{1};
data.drug_redfern = C{2};
data.drug_IC50_Na = C{3};
data.drug_IC50_Ca = C{4};
data.drug_IC50_HERG = C{5};
data.drug_low_ETPC = C{6};
data.drug_high_ETPC = C{7};
data.drug_in_redfern_figures = C{8};
data.grandi_measure = C{9};

% Convert IC50 values of -2 (no data) into very large IC50 values 
% (no/very small effect on current) and don't mess up classification.
indices = abs(data.drug_IC50_Na+2) < 1e-8;
data.drug_IC50_Na(indices) = 1e10;
indices = abs(data.drug_IC50_Ca+2) < 1e-8;
data.drug_IC50_Ca(indices) = 1e10;
indices = abs(data.drug_IC50_HERG+2) < 1e-8;
data.drug_IC50_HERG(indices) = 1e10;


% This is the concentrations that all drugs are run at (in nM).
data.drug_concentrations_text = {'0','1','10','100','1000','10000','100000','1e+06'};
data.drug_concentrations = [0;1;10;100;1e3;1e4;1e5;1e6];

data.drug_IC50_Na_over_HERG = data.drug_IC50_Na./data.drug_IC50_HERG;
data.drug_IC50_Ca_over_HERG = data.drug_IC50_Ca./data.drug_IC50_HERG;

for i=1:length(data.drug_redfern)
    if data.drug_redfern(i) == 1
        redfern_colour(i) = 'k';
    elseif data.drug_redfern(i) ==2
        redfern_colour(i) = 'r';
    elseif data.drug_redfern(i) == 3 
        redfern_colour(i) = 'm';
    elseif data.drug_redfern(i) == 4
        redfern_colour(i) = 'b';
    elseif data.drug_redfern(i) == 5 
        redfern_colour(i) = 'g';   
    else
        error('Unrecognised Redfern category')
    end
end
data.drug_redfern_colour = redfern_colour;

data.drug_indices = [0:length(data.drug_names)-1]';
