%compare amount of time between the CPU and GPU 
%for each step

clear, clc, close all;

A1 = gpuArray(rand(1200));
A2 = gpuArray(rand(1200));

B1 = rand(1200);
B2 = rand(1200);

num_imgs = 100000;

disp('background subtraction')
disp('gpu')
tic
for i = 1:num_imgs
    A3 = A1 - A2;
end
toc/num_imgs
a = toc/num_imgs;

disp('cpu')
tic
for i = 1:num_imgs
    B3 = B1 - B2;
end
toc/num_imgs
b = toc/num_imgs;

disp('-------------------------------------')
disp('median filter')
disp('gpu')
tic
for i = 1:100
    A4 = medfilt2(A3,[9 9]);
end
toc/100
a = a + toc/100;

disp('cpu')
tic
for i = 1:100
    B4 = medfilt2(B3,[10 10]);
end
toc/100
b = b + toc/100;
disp('-------------------------------------')
disp('find circles')
disp('gpu')
tic
for i = 1:100
    [A5, A6] = imfindcircles(A4,[10 20],'ObjectPolarity','bright',"Sensitivity",.95);
end
toc/100
a = a + toc/100;

disp('cpu')
tic
for i = 1:100
    [B5, B6] = imfindcircles(B4,[10 20],'ObjectPolarity','bright',"Sensitivity",.95);
end
toc/100
b = b + toc/100;
disp('-------------------------------------')
disp('time with GPU: ' + string(a) + 's per img')
disp('time with GPU: ' + string(a*100000/3600) + 'h per 100k img')
disp('-------------------------------------')
disp('time with CPU: ' + string(b) + 's per img')
disp('time with CPU: ' + string(b*100000/3600) + 'h per 100k img')
