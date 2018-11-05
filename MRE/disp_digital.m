%% display all slice number with digital number output
function slice_index = disp_digital (Nx, Ny, sliceNo)
clear slice_index
% slice_index = zeros(Nx,Ny,sliceNo);
for i = 1: sliceNo
    a=sevenSegmentDisplay(floor(i/10));
    a1=a.sevSegDispMatrix(:,:,2)./a.black;
    a=sevenSegmentDisplay(mod(i,10));
    a2=a.sevSegDispMatrix(:,:,2)./a.black;
    tmp = cat(2,a1,a2);
    if(round((Nx-size(tmp,1))/2)<0)
        tmp2=tmp;
        slice_index(:,:,i) = single(tmp2);
    else
        tmp2 = padarray(tmp,[round((Nx-size(tmp,1))/2),round((Ny-size(tmp,2))/2)],'both');
        slice_index(:,:,i) = single(tmp2(1:Nx,:));
    end
    
end