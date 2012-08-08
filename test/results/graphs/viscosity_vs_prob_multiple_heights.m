%function [ output_args ] = viscosity_vs_prob_multiple_heights( input_args )
%PERCENT_FIXED_VS_HEIGHT Summary of this function goes here
%   Detailed explanation goes here
[data, settings] = loadData();
num_heights = length(settings.heights);

% These are hard-coded into the simulations
crypt_length = 14.500;
box_proportion = 1.0/20.0;
data_to_analyse=[];
for v=1:length(settings.viscosities)
   % figure
    for m=1:length(settings.mutations)
        % Only load 'initial mutation height'
        dataset = data{m,1,v};
%         for h=1:length(settings.heights)
%             % For all heights and mutations (viscosity 2 = 1.0, control)
%             dataset = [dataset;data{m,h,v}];
%         end
        data_to_analyse{m,v} = dataset;
    end
end

clear dataset
b=1;
Success=[];
Total=[];
plotting={[],[],[],[]};
plotting_viscos={[],[],[],[]};
num_bucket={[],[],[],[]};
for v = 1:length(settings.viscosities)
    for m = 1:length(settings.mutations)
        dataset = data_to_analyse{m,v};
        if length(dataset) > 600
            dataset = dataset(1:600,:);
        end
        successes=0;
        failures=0;
        
        for i=1:length(dataset)
            
            if dataset(i,4)==1
                successes(b) = successes(b)+1;
            else
                failures(b) = failures(b) +1;
            end
        end
        
        assert(sum(successes)+sum(failures)==length(dataset));
        
        disp(['Dataset v = ' settings.viscosity_names{v} ', m = ' settings.mutation_names{m}, ' #sims = ' num2str(successes+failures)])
        
        % Condition on viscosities put in here to avoid stretching errorbar
        % plots below.
        if (successes+failures>0 && settings.viscosities(v)<=3.9)
            Success(m,v) = successes;
            Total(m,v) = successes+failures;
            Proportions(m,v) = 100.*successes./(successes+failures);
            plotting{m} = [plotting{m};Proportions(m,v)]
            plotting_viscos{m} = [plotting_viscos{m};settings.viscosities(v)]
            num_bucket{m} = [num_bucket{m};Total(m,v)];
        end
    end
end

% Order in increasing proliferation height
order = [2 4 3 1];
prolifs = [50 100 90 35];
for m=1:length(settings.mutations)
    figure(123)
    
    % Use this to print all 4 y_thr values.
    %subplot(1,4,order(m))
    
    % Use these lines to just print y = 0.5, 0.9, 1.0
    if (m==4)
        continue
    end
    subplot(1,3,order(m)-1)
    
    plot(plotting_viscos{m},plotting{m},'k--')
    hold on
    set(gca,'FontSize',14)

    ylim([0 50])
    xlabel('Mutant adhesion parameter \alpha','FontSize',16)
    ylabel('Domination probability (%)','FontSize',16)
    title(['y_{thr} = ' num2str(prolifs(m)) '%'],'FontSize',16)
    
    confidence_interval = 0.975; % 95% CI altogether for both sides
    % for standard normal (central limit thm)
    errors = norminv(confidence_interval,0,1) .* sqrt(((plotting{m}./100).*(1.0-plotting{m}./100))./num_bucket{m});
    
    errorbar(plotting_viscos{m},plotting{m},100*errors,'kx')
end


% %% Export this figure to the paper images directory.
% papersize = get(gcf, 'PaperSize');
% width = 30.0;         % Initialize a variable for width.
% height = 12.0;      % Initialize a variable for height.
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(gcf, 'PaperPosition', myfiguresize);

%print -depsc2 ../../../paper/images/proliferation_vs_viscosity



