function U = interp(u,X)
global nv h;
U = zeros(nv,3);
for k = 1:nv
  s = X(k,:)/h; 
  i = floor(s); % integer part of s
  r = s-i;      % decimal part of s
  i1 = ((i(1)-1):(i(1)+2))+1;    
  i2 = ((i(2)-1):(i(2)+2))+1;
  i3 = ((i(3)-1):(i(3)+2))+1;
  
  % w has 64 combinations of phi because each phi are 4 points supported.
  ww = phi1(r(1)).*phi2(r(2)).*phi3(r(3));
  U(k,1) = sum(sum(sum(ww.*u(i1,i2,i3,1))));
  U(k,2) = sum(sum(sum(ww.*u(i1,i2,i3,2))));
  U(k,3) = sum(sum(sum(ww.*u(i1,i2,i3,3))));
end