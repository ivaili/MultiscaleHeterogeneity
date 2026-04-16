function out = flattenMat(a,k)
if nargin<2
 k=1;
end
mask = triu(true(size(a)),k);
out = a(mask);
end
%k=0 includes diagonal; k=1 doesn't
