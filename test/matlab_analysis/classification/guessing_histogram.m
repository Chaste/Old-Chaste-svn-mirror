% What would the plot look like if we just guessed at random?!

num_samples = 1e6;

actual_group = randi(5,num_samples,1);
guess_group = randi(5,num_samples,1);

errors = abs(actual_group-guess_group);


figure(1)
hist(errors)
xlabel('Category Error')
ylabel('Frequency')

buckets = zeros(5,1);
for i=1:num_samples
    for j=1:5
        if errors(i)==j-1
            buckets(j) = buckets(j)+1;
        end
    end    
end

proportions = buckets./num_samples

mean_error = mean(errors)
std_error = std(errors)

ideal_proportions = [5/25 8/25 6/25 4/25 2/25]'
ideal_mean = 40/25
ideal_std = sqrt(100/25 - (40/25).^2)