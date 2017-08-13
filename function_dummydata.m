function [ ims imdepths ] = function_dummydata( Setup )
%%% This function for the purposes only of exemple simulation, provides a
%%% series of masks (positive) and depths for the desired format. 

z = [-300 : 2: 300] * 30e-4 ;   % Depth level requested in 3D region.
load('Letters_Data');
masks = masks(200:500,200:500,:);
n = 10;
% % Preprocess image
resizeX =450;
resizeY = 600;
cropxmin = 1; cropxmax=450; cropymin=1; cropymax=600;
shiftx = -50;
shifty = 0;
padx = 80;
pady = 80;
masks = imresize(masks, [resizeX, resizeY]);
masks = circshift(masks(cropxmin:cropxmax, cropymin:cropymax, :), [shiftx, shifty, 0]);
masks = padarray(masks, [padx, pady, 0]);
ims = imresize(permute(masks(:,:,1:n), [2,1,3]) ,[Setup.Nx, Setup.Ny]);
step = floor(50 / n) + 1;
indices = [21:step:71];
if numel(indices) < n
    indices(end+1) = indices(end) + step;
end
imdepths = z(indices);
z = z(indices);
nfocus = floor(numel(z)/2)+1;
focus = (z(end) + z(1))/2;
zdisplay = imdepths;
ims = max(ims,0);
if Setup.verbose == 1;
disp('Dummy data Ready !')
end
imdepths = imdepths-mean(imdepths);
imdepths = 6*imdepths;


ims = double(ims>0.3);
end

