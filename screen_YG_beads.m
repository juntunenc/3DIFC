clc, clear, close all;

% screen YG beads

img_folder = 'C:\Users\juntune3\Desktop\20221111YG_Beads\2';

bg_img = 'C:\Users\juntune3\Desktop\20221111YG_Beads\2\11-41-29.000-0002.tif';


img_mat = screen_imgs(img_folder,bg_img,false,.9);

for i = 1:round(size(img_mat,1)/100)
    start_img = ((i-1)*100)+1;
    end_img = i * 100;
    img_mat_save = img_mat(start_img:end_img,:,:);
    save(['C:\Users\juntune3\Desktop\20221111YG_Beads/img_mat'+string(i)+'.mat'],'img_mat_save');
end

