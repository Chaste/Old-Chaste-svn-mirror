%% Make a figure of how good predictions tend to be...
% Requires script.m to have been run for 1D and 2D classification.
close all
clc

load('compiled_results.mat');
original_data = compiled_results.measures;
original_names = compiled_results.measure_names';

chosen_marker_index = 316;

% For 2D classification problem stick to just examining multi/herg-only
% cases - not mixtures. HARDCODED VALUES (when extra measures are added
% these need changing).
D1_multi_channel_indexes = [1 2 6 7 9 10 12:386];
D1_herg_only_indexes = [3:5 8 11 387:length(compiled_results.measure_names)];

remove_cat_1 = false;
combine_cats_1_and_2 = true;

if remove_cat_1
    disp(['For classification involving categories ' num2str(remove_cat_1+1) '-5'])
    extra_string = '_no_cat_1';
elseif combine_cats_1_and_2
    disp(['For classification involving categories 1&2, 3, 4, 5'])
    extra_string =  'combine_12';
else
    extra_string = '';
end

for dimension_of_analysis = 1:2
    for stratification=5:-1:1
        figure
        if stratification==5
            load(['mean_error_in_category_' num2str(dimension_of_analysis) 'D_' extra_string '.mat'])
            load(['std_dev_error_in_category_' num2str(dimension_of_analysis) 'D_' extra_string '.mat'])
        else
            load(['mean_error_in_category_' num2str(dimension_of_analysis) 'D_strat_' num2str(stratification) '_' extra_string '.mat'])
            load(['std_dev_error_in_category_' num2str(dimension_of_analysis) 'D_strat_' num2str(stratification) '_' extra_string '.mat'])
        end
        
        if dimension_of_analysis == 1 % We want to compare with original 1D markers
            meanHerg(stratification) = mean_error_in_category(8);
            stdHerg(stratification) = std_dev_error_in_category(8);
            meanRedfern(stratification) = mean_error_in_category(11);
            stdRedfern(stratification) = std_dev_error_in_category(11);
            meanDecentMarker(stratification) = mean_error_in_category(chosen_marker_index);
            stdDecentMarker(stratification) = std_dev_error_in_category(chosen_marker_index);
        end
        
        fprintf('Mean Error: %g, std: %g, measure: %s \n',...
            meanHerg(stratification),stdHerg(stratification),original_names{8} )
        
        fprintf('Mean Error: %g, std: %g, measure: %s \n',...
            meanRedfern(stratification),stdRedfern(stratification),original_names{11} )
        
        % Convert 1 and 2D measures into a long list we can sort...
        counter = 0;
        herg_only_indexes = [];
        multi_channel_indexes = [];
        clear first_measure second_measure mean_error_in_cat std_dev_error_in_cat
        for i=1:size(mean_error_in_category,1)
            for j=1:size(mean_error_in_category,2)
                if mean_error_in_category(i,j)>0
                    counter = counter + 1;
                    first_measure(counter) = j;
                    second_measure(counter) = i;
                    mean_error_in_cat(counter) = mean_error_in_category(i,j);
                    std_dev_error_in_cat(counter) = std_dev_error_in_category(i,j);
                    if dimension_of_analysis==2
                        if (length(find(D1_herg_only_indexes==i))>0 && ...
                                length(find(D1_herg_only_indexes==j))>0)
                            % This is a completely "herg only measure"
                            herg_only_indexes = [herg_only_indexes counter];
                        else
                            multi_channel_indexes = [multi_channel_indexes counter];
                        end
                    end
                end
            end
        end
        
        if dimension_of_analysis==1
            multi_channel_indexes = D1_multi_channel_indexes;
            herg_only_indexes = D1_herg_only_indexes;
        end        
        
        % Remove repeated data points so that the eps files aren't
        % ma-hussive.
        [things_to_plotx things_to_ploty] = remove_repeated_points(mean_error_in_cat(multi_channel_indexes), std_dev_error_in_cat(multi_channel_indexes));
        scatter(things_to_plotx, things_to_ploty, 'ko')
        hold on
        [things_to_plotx things_to_ploty] = remove_repeated_points(mean_error_in_cat(herg_only_indexes), std_dev_error_in_cat(herg_only_indexes));
        scatter(things_to_plotx, things_to_ploty, 'k*')
        
        % Try to say which measures are best
        [sorted_mean_errors idx] = sort(mean_error_in_cat + 1e-7*std_dev_error_in_cat);% Cheeky way of sub sorting over std for same mean...
        fprintf('Top %iD Classifiers by mean error:\n',dimension_of_analysis)
        fprintf('mean & standard deviation & measure\n')
        for i=1:30
            fprintf('$%4.3f$ \t& \t $%4.3f$ \t & %s ',...
                mean_error_in_cat(idx(i)),std_dev_error_in_cat(idx(i)),original_names{first_measure(idx(i))} )
            if dimension_of_analysis==2
                fprintf(' and %s', original_names{second_measure(idx(i))})
            end
            % See which points on the graphs these correspond to
            if stratification==5
                % Plot the top 30 predictors in bold
                scatter(mean_error_in_cat(idx(i)),std_dev_error_in_cat(idx(i)),'ko','LineWidth',2.0)
                best_predictors_(i) = idx(i);
            else % See where these are in the stratified groups.
                scatter(mean_error_in_cat(best_predictors_(i)), std_dev_error_in_cat(best_predictors_(i)),'ko','LineWidth',2.0)
            end    
            hold on
            fprintf('\\\\ \n')
        end
        fprintf('\n')
               
        if remove_cat_1 || combine_cats_1_and_2
            expected_mean = 20/16;
            expected_std = sqrt(40/16 - (20/16).^2);
        else
            expected_mean = 40/25;
            expected_std = sqrt(100/25 - (40/25).^2);
        end
        
        if stratification < 5
            title_suffix = ['Stratification ' num2str(stratification)];
        else
            title_suffix = '';
        end
                
        title(title_suffix,'FontSize',16)
%         if remove_cat_1
%             title(['Categories 2-5' title_suffix],'FontSize',16)
%         elseif combine_cats_1_and_2
%             title(['Categories 1&2,3,4,5' title_suffix],'FontSize',16)
%         else
%             title(['Categories 1-5' title_suffix],'FontSize',16)
%         end
        

        %if dimension_of_analysis==1
            plot([expected_mean expected_mean],[0 3],'k-')
            plot([0 5],[expected_std expected_std],'k-')
            plot([meanRedfern(stratification) meanRedfern(stratification)],[0 3],'k--')
            plot([0 5],[stdRedfern(stratification) stdRedfern(stratification)],'k--')
            plot([meanDecentMarker(stratification) meanDecentMarker(stratification)],[0 3],'k:')
            plot([0 5],[stdDecentMarker(stratification) stdDecentMarker(stratification)],'k:')
            set(gca,'FontSize',14);
            xlim([0 2.5])
            ylim([0 1.5])
            xlabel('Mean category prediction error','FontSize',16)
            ylabel('Standard deviation of category prediction error','FontSize',16)
        %end
        %legend('1D classification', '2D classification', 'Expected values for uniform distribution','Location','SouthEast')
    end
end


