function [ ims, imdepths ] = function_dummydata( Setup )
% % [ ims, darkims, imdepths ] = function_dummydata( Setup )
%%% This function for the purposes only of exemple simulation, provides a
%%% series of target intensity (ims), dark voxels (darkims), and depths (imdepths) for the desired format. 

load('Dummy_data.mat');
imdepths =0.1*linspace(1,10,10);
masks = masks(220:380,325:485,:);
ims = imresize(masks, [Setup.Nx, Setup.Ny]);

if Setup.verbose == 1;
disp('Dummy data Ready !')
end

ims = double(ims>0.3);

v = squeeze(sum(sum(double(ims>0.1),1),2)); v = v/(max(v)); v = 1./v;
for i = 1:10
    pick = squeeze(ims(:,:,i)); 
    pick = imresize(pick, sqrt(v(i)));
    [cx cy] = size(pick);
    pick = pick(floor(cx/2-Setup.Nx/2+1):floor(cx/2+Setup.Nx/2),floor(cy/2-Setup.Ny/2+1):floor(cy/2+Setup.Ny/2));
        pause(0.01)
 ims(:,:,i) =pick;   
end
ims = double(ims>0.1);
v = squeeze(sum(sum(double(ims>0),1),2)); 


LN = 5;
ims = ims(:,:,1:LN);
imdepths = imdepths(1:LN);
imdepths =imdepths -mean(imdepths);

darkims = 40*ims;
ims(:,:,2) = 0;
ims(:,:,4) = 0;
darkims(:,:,1) = 0;
darkims(:,:,3) = 0;
darkims(:,:,5) = 0;

%darkims = darkims-0.01;


end

