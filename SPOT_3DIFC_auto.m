%% Code to process the 3DIFC data 
% Shreyas start here
%%
clear
close all
tic

global res kres padd0
global row0 col0 win0 Nx Ny
global Pmla F1 F2 Ftube Mcond n0 k0n


%% Load the FTS-processed data

datapath = strcat('C:\Users\juntune3\Desktop\Projects\3DIFC\10.18.22\750us_1ms_6');

img = double(imread([datapath '\13-29-51.000-1.tif']));
bg = double(imread([datapath '\13-29-51.000-2.tif']));

img = img - bg;

img(img<0) = 0;
figure(1), imagesc(img); axis image;

filt_img = medfilt2(img,[10 10]);

figure(99), imagesc(filt_img);
%% Get image centered

[centers radii] = imfindcircles(filt_img,[10 14], Sensitivity=.95);

h = viscircles(centers,radii);

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

figure(100), imagesc(lines) %this is for visualizing where it thinks the beads are

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
figure(101),imagesc(cen_img)


V = [round(size(img,1)/2)-y_cen round(size(img,2)/2)-x_cen];
img = circshift(img, V);
figure(102),imagesc(img)

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


Showgrid(img);
toc

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


figure(51), imagesc(proj_img), axis image, axis off; colormap gray;

figure(52), imagesc(log10(abs(squeeze(F_tomo3(:,:,padd0/2+1)))));
axis image; colormap jet; colorbar; set(gca,'fontsize',13); axis off;

figure(53), imagesc(log10(abs(squeeze(F_tomo3(:,padd0/2+1,:)))));
axis image; colormap jet; colorbar; set(gca,'fontsize',13); axis off;

% Regularization
% F_tomo3 = Regularization(F_tomo3);
% F_tomo3 = HandleSingularity(F_tomo3);
tomo_f = ifftshift(ifftn(fftshift(F_tomo3)))*(kres*padd0)^3;
%

I3d = real(tomo_f);

%% Show the horizontal cross-sections at different heights
kcen = padd0/2+1;
for kk2 = kcen-50:kcen+50
    %     kk2 = padd0/2 + 1;
    figure(61), imagesc(I3d(:,:,kk2)');
    axis image; title(['kk = ' num2str(kk2)]); colormap gray;
    pause(0.1);
end

%% Show the center horizontal and vertical cross-sections
% I3d = I3d/max(I3d(:));

kk2 = kcen;
figure(62), imagesc(I3d(:,:,kk2));
axis image; colormap jet; axis off; colorbar;
set(gca,'fontsize',18);

ii2 = 86;
xtick = ((1:length(I3d))-ii2)*res*1e6;
figure(63), plot(xtick,I3d(ii2,:,kk2),'k-+'); grid;
xlabel('Coordinate along the dotted line (µm)');
xlim([-10 10]); set(gca,'xtick',[-10:2:10]); set(gca,'fontsize',16);

figure(64), imagesc(flipud(squeeze(I3d(ii2,:,:)))');
axis image; colormap jet; axis off; colorbar; colormap gray;
set(gca,'fontsize',18);

%% 3-D rendering
padd1 = padd0;

%I3d = I3d - min(I3d,[],'all');
%I3d = I3d / max(I3d,[],'all');

temp1 = sum(I3d,3);
temp1 = temp1 - min(temp1,[],'all');
temp1 = temp1/max(temp1,[],'all');
figure(71), imagesc(temp1)

mask1 = zeros(size(temp1));

mask1(temp1>.4) = 1; %mask to change based on figure 71
figure(72), imagesc(mask1)
mask1 = imdilate(mask1,strel('disk',3));
mask1 = imerode(mask1,strel('disk',5));
mask1 = imdilate(mask1,strel('disk',7));
figure(73), imagesc(mask1)

%temp2 = squeeze(sum(I3d,1));
%figure(74), imagesc(temp2)
% mask2 = zeros(size(temp2));
%
% mask2(temp2>20) = 1;
% mask2 = imerode(mask2,strel('disk',5));
% mask2 = imdilate(mask2,strel('disk',5));
%
% mask = repmat(mask1,[1 1 padd1]).*permute(repmat(mask2,[1 1 padd1]),[3 1 2]);
I3d_masked = I3d.*mask1;
I3d_masked = circshift(I3d_masked,[0 0 0]);

xtick = ((0:padd0-1)-padd0/2)*res/1e-6;
[X,Y,Z] = meshgrid(xtick,xtick,xtick);

xmin = -10; xmax = 10;
zmin = -10; zmax = 10;

xslice = [];
yslice = [];
zslice = [zmin:1:zmax];

figure(2000);
h = slice(X,Y,Z,I3d_masked,xslice,yslice,zslice);

xlim([xmin xmax]); ylim([xmin xmax]); zlim([zmin zmax]);
xticks([xmin:5:xmax]); yticks([xmin:5:xmax]); zticks([zmin:2:zmax]);

xlabel(['X (' char(181) 'm)']); ylabel(['Y (' char(181) 'm)']); zlabel(['Z (' char(181) 'm)']);

set(gca,'fontsize',18);

colormap(flipud(gray));
%colormap((gray));

set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp');

alpha('color');

%alphamap('default');
alphamap('rampup');
alphamap('decrease',.4);
view([1 -1 .4]);

%OptionZ.FrameRate=30;OptionZ.Duration=10;OptionZ.Periodic=true;
%CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'WellMadeVid',OptionZ)
%%
%view([1 -1 0.4]);
%
% test_slice = I3d(:,:,90);
% figure(),imagesc(test_slice)
%
% mask = zeros(size(I3d,1),size(I3d,2),size(I3d,3));
% for i = 1:size(I3d,3)
%     i
%     temp1 = I3d(:,:,i);
%     temp1 = temp1 - min(temp1,[],'all');
%     temp1 = temp1/max(temp1,[],'all');
%     %figure(71), imagesc(temp1)
%     %pause()
%     temp_mask = zeros(size(temp1));
%     temp_mask(temp1<.01) = 1;
%     mask(:,:,i) = temp_mask;
%
%     %figure(72),imagesc(temp_mask)
%     %pause()
% end
%
% I3d_masked = I3d.*mask;
%
% xtick = ((0:padd0-1)-padd0/2)*res/1e-6;
% [X,Y,Z] = meshgrid(xtick,xtick,xtick);
%
% xmin = -20; xmax = 20;
% zmin = -20; zmax = 20;
%
% xslice = [];
% yslice = [];
% zslice = [zmin:5:zmax];
%
% figure(2001);
% h = slice(X,Y,Z,I3d_masked,xslice,yslice,zslice);
%
% xlim([xmin xmax]); ylim([xmin xmax]); zlim([zmin zmax]);
% xticks([xmin:5:xmax]); yticks([xmin:5:xmax]); zticks([zmin:2:zmax]);
%
% xlabel(['X (' char(181) 'm)']); ylabel(['Y (' char(181) 'm)']); zlabel(['Z (' char(181) 'm)']);
%
% set(gca,'fontsize',18);
%
% colormap(flipud(gray));
% %colormap((gray));
%
% set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp');
%
% alpha('color');
%
% %alphamap('default');
% alphamap('rampup');
% alphamap('increase',.5);
%
% %view([1 -1 0.4]);



% %%
% clear multi
%
% kk2 = 130;
% zslice = [-6:3:6];
% zz = kk2 + round(zslice/(res/1e-6));
%
% multi = I3d(:,:,zz(1));
% for ctr = 2:length(zz)
%     tmp = I3d(:,:,zz(ctr));
%     multi = cat(3,multi,tmp);
% end
%
% figure(53), montage(multi,'Size',[1 length(zz)],'DisplayRange',[0 50]);
% colormap(gray); colorbar; set(gca,'fontsize',18);
%
% %%
% kk2 = kcen;
% tmp = I3d(:,:,kk2);
%
% figure(81), imagesc(tmp), axis image; axis off;
%
% tmp1 = tmp>10;
% se = strel('disk',5);
% tmp2 = imerode(tmp1,se);
% tmp3 = imdilate(tmp2,se);
%
% figure(83), imagesc(tmp3), axis image; axis off;
%
% tmp4 = tmp(tmp3>0);
% mean(tmp4(:))


