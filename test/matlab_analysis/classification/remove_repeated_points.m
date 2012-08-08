function [things_to_plotx, things_to_ploty] = remove_repeated_points(input_vec1, input_vec2)
% Make sure we only put points on the scatter plot that haven't already
% been plotted (to stop us having eps files with a million points on them.)
things_to_plotx = input_vec1(1);
things_to_ploty = input_vec2(1);
N = length(input_vec1);
for i=2:N
    %fprintf('%i : Mean : %f, std: %f\n',i,mean_error(i),std_error(i))    
    x_match = find(abs(things_to_plotx(:) - input_vec1(i)) < 1e-6);
    if ~isempty(x_match)
        y_match = find(abs(things_to_ploty(x_match) - input_vec2(i)) < 1e-6);
        got_already = ~isempty(y_match);
    else
        got_already = false;
    end    
    if ~got_already
        things_to_plotx = [things_to_plotx; input_vec1(i)];
        things_to_ploty = [things_to_ploty; input_vec2(i)];
    end    
end