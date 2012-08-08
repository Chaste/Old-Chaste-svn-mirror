function h = redfern_gscatter_wrapper(x,y,g,marker_size,override_marker,group_to_redfern_mapping)
%GSCATTER   Gary's wrapped version of gscatter
% Has the distinct advantage of giving the same groups
% the same freaking colours and symbols each time...

if isempty(x) % Lets us use either training_points{cell} or x,y.
    xd=[];
    yd=[];
    gd=[];
    for i=1:length(y)
        temp = y{i};
        xd = [xd; temp(:,1)];
        yd = [yd; temp(:,2)];
        gd = [gd; ones(size(temp,1),1)*i];
    end
    x = xd;
    y = yd;
    g = gd;
end

if nargin>=6
    g = group_to_redfern_mapping(g);
end

background_shading=false;
% override symbols for the meshgrid shading of the entire plot.
if nargin>=5 && ~isempty(override_marker)
    background_shading=true;
end

colour = [];
symbol = [];

test_hold = ishold();

for i=1:5
    idx = find(g==i);
    if i==1
        colour = 'k'; % Add black
        symbol = 'v';
    elseif i==2
        colour = 'r'; % Add red
        symbol = 's';
    elseif i==3
        colour = 'm'; % Add magenta
        symbol = 'x';
    elseif i==4
        colour = 'b'; % Add blue
        symbol = '+';
    elseif i==5
        colour = 'g'; % Add green
        symbol = '*';
    end
    if background_shading
        symbol = override_marker;
    end
    
    h = scatter(x(idx),y(idx),[colour symbol],'SizeData',marker_size,'LineWidth',2);
    hold on
end

if ~test_hold
    hold off % Only switch the hold off if it was off before the call to this method.
end

