%%
function Showgrid(img)

global row0 col0 win0 Nx Ny

[xsize,ysize] = size(img);

xi = row0:win0:row0+win0*Nx;
yi = col0:win0:col0+win0*Ny;

h = figure(1);
if ishandle(h)
    close(h)
    h = figure(1);
end

imagesc(img,[min(img(:)) max(img(:))*0.5]), axis image, colormap gray; hold on;
movegui(h,'northwest');

for m = 1:Nx+1
    line([1 ysize],[xi(m) xi(m)],'color','white');
end

for n = 1:Ny+1
    line([yi(n) yi(n)],[1 xsize],'color','white');
end

axis([1 ysize 1 xsize])