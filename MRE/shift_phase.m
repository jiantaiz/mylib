function M2 = shift_phase(M, ang)
Re = real(M);
Im = imag(M);

if numel(ang)<2
    ang(2) = ang(1);
end

M2 = angle(exp(1i*(Re + ang(1)))) + 1i.*angle(exp(1i*(Im + ang(2)))) ;
