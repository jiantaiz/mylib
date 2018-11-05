function r = bsx(a,op,b)
sz_a = size(a);
sz_b = size(b);
if length(sz_b)<length(sz_a)
    sz_b(length(sz_b)+1:length(sz_a))=1;
end

rep = sz_a ./sz_b;
b = repmat(b,rep);
switch op
    case '+'
        r = a+b;
    case '-'
        r = a-b;
    case {'*', '.*'}
        r = a.*b;
    case {'/', './'}
        r = a./b;
end

