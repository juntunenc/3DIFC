function vol = getVol(img,vol_thresh)

%% Code to process the 3DIFC data 
% Shreyas start here
%%

global res kres padd0
global row0 col0 win0 Nx Ny
global Pmla F1 F2 Ftube Mcond n0 k0n


%% Load the FTS-processed data

img(img<0) = 0;
%figure(1), imagesc(img); axis image;

filt_img = medfilt2(img,[10 10]);

%figure(99), imagesc(filt_img);
%% Get image centered

[centers radii] = imfindcircles(filt_img,[10 14], Sensitivity=.95);

%h = viscircles(centers,radii);

%find the average xval of centers within 10 pixels of eachothers x
%coordinate
xs = centers(:,1);
ys = centers(:,2);

[xs_sorted, I_xs_sorted] = sort(xs);
[ys_sorted, I_ys_sorted] = sort(ys);

x_group_num = 1;

for i = 1:length(xs)
    
    if i == 1
        xs_sorted(1,2) = 1;

    elseif xs_sorted(i,1) < (xs_sorted(i-1,1) + 20)
        xs_sorted(i,2) = x_group_num;

    else
        x_group_num = x_group_num + 1;       
        xs_sorted(i,2) = x_group_num;

    end
end

%average the values of each group
for i = 1:x_group_num
    x_means(i) = round(mean(xs_sorted(xs_sorted(:,2) == i)));
    
end

y_group_num = 1;

for i = 1:length(ys)
    
    if i == 1
        
        ys_sorted(1,2) = 1;

    elseif ys_sorted(i,1) < (ys_sorted(i-1,1) + 20)
        ys_sorted(i,2) = y_group_num;

    else
        y_group_num = y_group_num + 1;       
        ys_sorted(i,2) = y_group_num;


    end
end

for i = 1:y_group_num
    y_means(i) = round(mean(ys_sorted(ys_sorted(:,2) == i)));
end

lines = img;
for i = 1:x_group_num
    lines(:,(x_means(i)-2):(x_means(i)+2)) = max(max(img));
end

for i = 1:y_group_num
    lines((y_means(i)-2):(y_means(i)+2),:) = max(max(img));
end

%figure(100), imagesc(lines) %this is for visualizing where it thinks the beads are

%find the center of the sample

if rem(length(x_means),2) == 0 %if even take average of middle two values
    x_cen = round((x_means(length(x_means)/2) + x_means(length(x_means)/2 + 1)) / 2);
else %else take middle value
    x_cen = median(x_means);
end
    
if rem(length(y_means),2) == 0 %if even take average of middle two values
    y_cen = round((y_means(length(y_means)/2) + y_means(length(y_means)/2 + 1)) / 2);
else %else take middle value
    y_cen = median(y_means);
end

%x_cen = round((x_means(1)+x_means(length(x_means)))/2);
%y_cen = round((y_means(1)+y_means(length(y_means)))/2);

cen_img = img;
cen_img(:,(x_cen-2):(x_cen+2)) = max(max(img));
cen_img((y_cen-2):(y_cen+2),:) = max(max(img));
%figure(101),imagesc(cen_img)


V = [round(size(img,1)/2)-y_cen round(size(img,2)/2)-x_cen];
img = circshift(img, V);
%figure(102),imagesc(img)

%% get window size

for i = 1:(length(x_means)-1)
    x_diffs(i) = x_means(i+1)-x_means(i);
end

for i = 1:(length(y_means)-1)
    y_diffs(i) = y_means(i+1)-y_means(i);
end

win_0_val = round((mean([mean(x_diffs) mean(y_diffs)]))/2)*2;
row_0_val = (y_means(1) - round(win_0_val/2)) + V(1);
col_0_val = (x_means(1) - round(win_0_val/2)) + V(2);

%Shreyas end here
%return win_0_val row_0_val col_0_val y_group_num col_0_val

%% Set the parameters to select projection images
% win0: the size of each projection images
% row0, col0: row and column index of the top, left corner
% Nx, Ny: the number of projection images along each dimension

win0 = win_0_val; % even number
row0 = row_0_val;
col0 = col_0_val;

Nx = y_group_num;
Ny = x_group_num;

padd0 = win0;


%Showgrid(img);
%toc

%% Set the parameters for SPOT processing (from the optical design)
Pmla = 500e-6; % MLA pitch (Edmund 64479)
Ftube = 200e-3; % Reference focal length (Nikon: 200e-3; Olympus: 180e-3; Zeiss: 165e-3)
Mcond = 60; % Objective lens magnification

Fobj = Ftube/Mcond;
F1 = 100e-3; % Focal length of Lens 1

F2 = 35e-3; % Focal length of Lens 2
Fmla = 13.8e-3; % Focal length of the micro-lens array (MLA)

Mag1 = (F1/Fobj)*(Fmla/F2);
% Mag2 = 10/200*50;
%Mag2 = 100e-3/40e-3*(75e-3/50e-3); % Magnification of the relay lenses (after the MLA)
Mag2 = 35e-3/20e-3;
Mag = Mag1*Mag2;

NAmla = Pmla/2/Fmla*Mag1;

pxl = 6.5e-6; % Camera pixel size
res = pxl/Mag;
xtick = ((0:padd0-1)-padd0/2)*res;

n0 = 1.337; % Refractive index of the medium
lambda = 500e-9;
k0 = 1/lambda;
k0n = 1/lambda*n0;

% Location of the optical axis in the raw image
% ii0 = 3.5, jj0 = 3.5 means the optical axis is located between the 3rd
% and 4th images along each axis
% rad0: the radius of the aperture (used to exclude the outermost images)
ii0 = y_group_num/2; jj0 = x_group_num/2; rad0 = 3;

%% SPOT processing
% This part requires understanding of the tomography reconstruction algorithm
kres = 1/(res*padd0);
ktick = ((0:padd0-1)-padd0/2)*kres;

[yyy0,xxx0] = meshgrid(1:padd0,1:padd0);

F_tomo3 = zeros(padd0,padd0,padd0);
N_rep3 = zeros(padd0,padd0,padd0);

proj_img = zeros(padd0,padd0);

flag = 0;

%iilist = circshift([1:Nx],-(ii0-1));
%jjlist = circshift([1:Ny],-(jj0-1));
iilist = circshift([1:Nx],-2);
jjlist = circshift([1:Ny],-2);

for ctr1 = 1:length(iilist)
    ii = iilist(ctr1);
    for ctr2 = 1:length(jjlist)
        jj = jjlist(ctr2);
        if (ii-ii0)^2+(jj-jj0)^2 <= rad0^2
            ij = (ii-1)*Ny + jj;
            row = row0 + (ii-1)*win0 + (win0-padd0)/2;
            col = col0 + (jj-1)*win0 + (win0-padd0)/2;

            % Align the projection images
            img_ij = img(row:row+padd0-1,col:col+padd0-1);

            sig = 1.85e-6/2.355/res;
            psf = fspecial('gaussian',30,sig);
            img_ij = deconvlucy(img_ij,psf);
            %             figure(72), imagesc(img_ij), axis image, colormap gray;

            if ~flag
                fixed = img_ij;
                flag = 1;
            else
                [optimizer, metric] = imregconfig('multimodal');
                img_ij = imregister(img_ij,fixed,'translation',optimizer,metric);
            end

            if mod(win0,2)
                img_ij = padarray(img_ij,[1 1],'replicate','post');
            end

            proj_img = proj_img + img_ij;

            % Calculation of the viewing angle (i.e., projection angle)
            % thx0, thy0: angles at the virtual plane wrt x and y axes
            % thx, thy: angles at the sample plane wrt x and y axes

            xxx = (ii-ii0)*Pmla*(F1/F2);
            yyy = (jj-jj0)*Pmla*(F1/F2);

            thx0 = atan(xxx/Ftube);
            thy0 = atan(yyy/Ftube);

            thx = asin(sin(thx0)*Mcond/n0); % Sine condition
            thy = asin(sin(thy0)*Mcond/n0);
            % Change starts from here
            [img_ij_proj,S] = Projection(img_ij,thx,thy);

            %             figure(61), subplot(1,2,1); imagesc(img_ij); axis image;
            %             subplot(1,2,2); imagesc(img_ij_proj); axis image;
            %             title([num2str(ii) ' ' num2str(jj)]);
            %             pause(0.5);

            F_img_ij = fftshift(fft2(ifftshift(img_ij_proj)))*res^2;

            Up = ktick(xxx0);
            Vp = ktick(yyy0);

            for pp = 1:padd0
                for qq = 1:padd0
                    UVW = S\[Up(pp,qq);Vp(pp,qq);0];
                    U = round(UVW(1)/kres) + padd0/2 + 1;
                    V = round(UVW(2)/kres) + padd0/2 + 1;
                    W = round(UVW(3)/kres) + padd0/2 + 1;
                    if U>1 && U<padd0 && V>1 && V<padd0 && W>1 && W<padd0
                        F_tomo3(U,V,W) = F_img_ij(pp,qq);
                    end
                end
            end
            % Change ends here
        end
    end
end



% Regularization
% F_tomo3 = Regularization(F_tomo3);
% F_tomo3 = HandleSingularity(F_tomo3);
tomo_f = ifftshift(ifftn(fftshift(F_tomo3)))*(kres*padd0)^3;
%

I3d = real(tomo_f);

%% calculate the volume
vox_size = res*10^6;

%
vol_norm = I3d - min(I3d,[],'all');
vol_norm = vol_norm / max(vol_norm,[],'all');

vol_mask = zeros(padd0,padd0,padd0);

vol_mask(vol_norm>vol_thresh)= 1;

%vol_mask(vol>.69)= 1; %sample 4 for 1um beads to match
vol = sum(vol_mask,"all")*vox_size^3; %um^3
%vol_sphere = 4/3*3.14*.5^3 %um^3


end