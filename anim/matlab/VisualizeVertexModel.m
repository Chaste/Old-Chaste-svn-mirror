function VisualizeVertexModel(filename)
% Function to visualize vertex cell model 
%
% Reads in <filename>.viznodes and <filename>.vizcells
%
% <filename>.viznodes has form:
% (time value 1) (x_1 y_1) (x_2 y_2) ... (x_m y_m)
% (time value 2) (x_1 y_1) (x_2 y_2) ... (x_n y_n)
% .
% .
% N.B. m not necessarily equal to n 
%
% <filename>.vizcells has form:
% (time value 1) (total number of nodes in cell 1) (node indices for cell 1) (total number of nodes in cell 2) (node indices for cell 2)....
%
% This requires LoadNonConstantLengthData('filename'), which is in the Chaste anim/matlab folder


nodesfile = [filename, '.viznodes'];
cellsfile = [filename, '.vizcells'];

nodedata = LoadNonConstantLengthData(nodesfile);
cellsdata = LoadNonConstantLengthData(cellsfile);

if size(nodedata) ~= size(cellsdata)

    error('Number of times does not match.')

end

% Plotting the nodes and elements

numtimes = length(nodedata);

for i = 1:numtimes

    time = nodedata{i}(1)     % Gives first element of ith row

    %numnodes = (length(nodes{i})-1)/2);

    xvals = nodedata{i}(2:2:end-1);
    yvals = nodedata{i}(3:2:end);

    hold off
    plot(xvals,yvals,'*')
    axis([0.9*min(xvals) 1.1*max(xvals) 0.9*min(yvals) 1.1*max(yvals)])
    title( ['Current time = ',num2str(time)] );
    hold on

    current_index = 2;          % Initialising

    while current_index < length(cellsdata{i})

        num_elem_nodes = cellsdata{i}(current_index);

        elem_nodes = cellsdata{i}(current_index+1:current_index+num_elem_nodes);     % Nodes for that cell       

        for j = 1:num_elem_nodes
            
            node1 = elem_nodes(j);
            node2 = elem_nodes(mod(j,num_elem_nodes)+1);
            xval1 = xvals(node1);
            xval2 = xvals(node2);
            yval1 = yvals(node1);
            yval2 = yvals(node2);
       
            plot([xval1, xval2], [yval1, yval2])

        end

        current_index = current_index + num_elem_nodes +1;           

    end

    pause(0.5);

end
