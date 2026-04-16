function [mat3d] = recon3dMat(flatMat,nodes,k)

if nargin<3
 k=1;
end
    
for i = 1:size(flatMat,2)
    mat3d(:,:,i) = reconMat(flatMat(:,i),nodes,k);
end
end

