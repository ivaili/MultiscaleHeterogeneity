function [flatmat] = flatten3dMat(matrix3d,k,csv_name)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin<2
 k=1;
end
for i = 1:size(matrix3d,3)
  a = matrix3d(:,:,i);
  flatmat(:,i)=flattenMat(a,k);

    
end
if nargin == 3
   dlmwrite(csv_name,flatmat,'delimiter',',') 
end

