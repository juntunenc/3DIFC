%%
function F_tomo3 = Regularization(F_tomo3)

global res kres padd0

F_tomo3 = HandleSingularity(F_tomo3);
F_mask = zeros(size(F_tomo3));
ind = F_tomo3~=0;
F_mask(ind) = 1;

f1 = ifftshift(ifftn(fftshift(F_tomo3)))*(kres*padd0)^3;

for ctr = 1:100
    F1 = fftshift(fftn(ifftshift(f1)))*res^3;
    
    f1 = real(f1);
    ind = f1<0;
    f1(ind) = 0;
        
    F1_new = fftshift(fftn(ifftshift(f1)))*res^3;
    
    F1p = F1.*F_mask + F1_new.*(1-F_mask);
    f1p = ifftshift(ifftn(fftshift(F1p)))*(kres*padd0)^3;
    f1 = f1p;
end

F_tomo3 = F1p;