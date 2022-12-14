function img_mat = screen_imgs_gpu(img_folder,bg_img,check_mode,sens)
tic
i = 0;

bg = gpuArray(double(imread([bg_img])));

tif_files = dir(fullfile(img_folder,'*.tif'));


for k = 1:length(tif_files)
    clc
    fprintf(1, 'file %s of %s \n', string(k),string(length(tif_files)));
    tot_time = toc;
    avg_time = tot_time/k;
    fprintf(1,'total elapsed time: %s s\n',string(tot_time))
    fprintf(1,'avg time per sample: %s s\n',string(avg_time))
    time_left = avg_time * (length(tif_files)-k) / 60;
    fprintf(1,'estimated: %s min remaining\n',string(time_left))

    baseFileName = tif_files(k).name;
    fullFileName = fullfile(img_folder, baseFileName);

    img = gpuArray(double(imread([fullFileName])));
    img = img - bg;

    filt_img = gpuArray(medfilt2(img,[9 9]));

    [centers, radii] = imfindcircles(filt_img,[10 20],'ObjectPolarity','bright',"Sensitivity",sens);

    if check_mode == true

        figure(1),imagesc(filt_img)
        viscircles(centers, radii,'Color','r');
        pause
    end


    if centers >= 20
        i = i + 1;
        img_mat(i,:,:) = img;

    end
    


end