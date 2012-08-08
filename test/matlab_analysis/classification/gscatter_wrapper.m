function h = gscatter_wrapper(x,y,g,varargin)
%GSCATTER   Gary's wrapped version of gscatter
% Has the distinct advantage of giving the same groups 
% the same freaking colours and symbols each time...

[g,sorted_idx] = sort(g);
x = x(sorted_idx);
y = y(sorted_idx);

h = gscatter(x,y,g,varargin{:});
