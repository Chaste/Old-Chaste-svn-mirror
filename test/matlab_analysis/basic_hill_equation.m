function percent_open = basic_hill_equation(conc, IC50)
hill = 1.0;

if IC50 < 0
    percent_open = 100.0;
else
    percent_open =  100.0/(1.0 + (conc/IC50)^hill);
end
