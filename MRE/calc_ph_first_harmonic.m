function out =  calc_ph_first_harmonic (cimg1,cimg2,tdim)
% CALC_PH_FIRST_HARMONIC calculates the first harmonic of a four phase offsets dataset by directly
% taking the 4-point FFT expression (x1+x3)+1i(x2+x4)
% out = CALC_PH_FIRST_HARMONIC(cimg1,[],tdim) returns the first harmonic of cimg1.
% out = CALC_PH_FIRST_HARMONIC(cimg1,cimg2,tdim) returns the first harmonic of cimg1.*conj(cimg2).


%
idx = repmat({':'},[1,numel(size(cimg1))]);
idx{tdim} = 1;
x1 = cimg1(idx{:});
idx{tdim} = 2;
x2 = cimg1(idx{:});
idx{tdim} = 3;
x3 = cimg1(idx{:});
idx{tdim} = 4;
x4 = cimg1(idx{:});

% x1 = cimg1(:,:,:,1,:);
% x2 = cimg1(:,:,:,2,:);
% x3 = cimg1(:,:,:,3,:);
% x4 = cimg1(:,:,:,4,:);

if isempty(cimg2)
    %     DC = exp(1i.*angle(x1.*x2.*x3.*x4)./4);
    %
    %     out = angle(x1.*conj(DC)) - angle(x3)
    
    out = angle(x1 .* conj(x3)) + angle(x2 .* conj(x4)).*1i;
else
    
        idx = repmat({':'},[1,numel(size(cimg2))]);
        idx{tdim} = 1;
        y1 = cimg2(idx{:});
        idx{tdim} = 2;
        y2 = cimg2(idx{:});
        idx{tdim} = 3;
        y3 = cimg2(idx{:});
        idx{tdim} = 4;
        y4 = cimg2(idx{:});
    %
%     y1 = cimg2(:,:,:,1,:);
%     y2 = cimg2(:,:,:,2,:);
%     y3 = cimg2(:,:,:,3,:);
%     y4 = cimg2(:,:,:,4,:);
%     
    
    out = angle( x1.*conj(y1) .* conj(x3).*y3) +  angle(x2 .* conj(y2) .* conj(x4) .*y4) .*1i;
end