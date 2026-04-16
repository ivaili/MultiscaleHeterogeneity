%reconstruct a connectivity matrix from a vector

function matrix = reconMat(vector,nNodes,k)
if nargin<3
 k=1;
end
matrix = zeros(nNodes);
upperInds = find(triu(ones(nNodes),k));
matrix(upperInds) = vector;
matrix = matrix+matrix';

if k == 0
matrix(logical(eye(size(matrix)))) = matrix(logical(eye(size(matrix))))/2;
else
    matrix(logical(eye(size(matrix))))=1;
end
end
%k=0 includes diagonal; k=1 doesn't
