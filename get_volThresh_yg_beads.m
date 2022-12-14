clear, clc, close all;

%load('C:/Users/juntune3/Desktop/20221024_YG_Beads/img_mat.mat');

num_files = 9;
tic
j = 0;

for k = 1:num_files
    
    load('C:\Users\juntune3\Desktop\20221111YG_Beads/img_mat' + string(k) + '.mat');
    img_mat = img_mat_save;
    clear img_mat_save
    for i = 1:size(img_mat,1)

        clc
        fprintf(1, 'file number: %s image %s of %s \n',string(k), string(i),string(size(img_mat,1)));
        tot_time = toc;
        avg_time = tot_time/(((k-1)*100)+i);
        fprintf(1,'total elapsed time: %s s\n',string(tot_time))
        fprintf(1,'avg time per sample: %s s\n',string(avg_time))
        time_left = avg_time * ((size(img_mat,1)-i)+(num_files-k)*100) / 60;
        fprintf(1,'estimated time for this file: %s min remaining\n',string(time_left))

        img = double(squeeze(img_mat(i,:,:)));

        try
            j = j + 1;
            [volThresh(j),vol(j)] = getVolThresh(img);
            disp(['Sample ' + string(i) + ': Threshold: ' + string(volThresh(j)) + ' Volume: ' + string(vol(j)) + ' um^3'])
        catch
            disp(['skipping sample: ' + string(i)])
        end

    end
end

toc

%volume of a 10um diameter sphere:
vol_sphere = 4/3*3.14*5^3; %um^3
disp(['Ground Truth: ' + string(vol_sphere) + ' um^3'])

figure(1),histogram(volThresh,30)
xlabel('Threshold Value')
ylabel('Occurrences')
figure(2),histogram(vol)
