echo 'apps
dealii
global
heart
linalg
mesh
ode' | while read dir; do
  find $dir -name "*.?pp" | while read file; do
    cat COPYING/file-template.txt $file > tmp
    mv tmp $file
  done
done
