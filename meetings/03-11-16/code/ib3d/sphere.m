%sphere.m: main script to generate a sphere with triangular mesh
nref = 4;
nvmax = 12*4^nref;
ntmax = 20*4^nref;
X = zeros(nvmax,3);
v = zeros(ntmax,3);
dodec
for nr=1:nref
  refine
end
X = X(1:nv,:);
