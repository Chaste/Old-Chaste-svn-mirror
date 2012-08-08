function new_points = transform_coordinates(old_points, transform_matrix)

assert(size(old_points,2)==size(transform_matrix,1))

for i=1:size(old_points,1)
    for j=1:size(old_points,2)
        new_points(i,j) = transform_matrix(:,j)'*old_points(i,:)';
    end
end