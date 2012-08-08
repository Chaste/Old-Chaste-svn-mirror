function APD = CalculateAPD( time, voltage, percentage )
%CALCULATEAPD Summary of this function goes here
%   Detailed explanation goes here

resting_voltage = voltage(1);
final_voltage = voltage(end);
voltage_range = max(voltage) - min(voltage);
threshold = resting_voltage+(max(voltage)-resting_voltage)*(100.0-percentage)/100.0;

if (abs(resting_voltage - final_voltage)>0.1*voltage_range)
    warning('Cell not at rest at beginning or end of trace');
    APD = -1;
    return;
end

under = true;
for i=1:length(voltage)
    if under==true
       if voltage(i) > threshold
           start_time = interpolate(time, voltage, threshold, i);
           under = false;
       end
    else
        if voltage(i) < threshold
            end_time = interpolate(time, voltage, threshold, i);
            APD = end_time - start_time;
            break
        end
    end
    if i==length(voltage)
        % No APD detected
        warning('No AP detected in this trace');
        APD = -1;
    end
end

function t = interpolate(time, voltage, threshold, i)
t = time(i-1) + ((threshold-voltage(i-1))/...
    (voltage(i)-voltage(i-1)))*(time(i)-time(i-1));