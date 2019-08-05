function B = imrotate3D(A,angle,rotax,varargin)
if nargin<3
    rotax = 'z';
end

switch rotax
    case 'x'
        A = permute(A,[1 3 2]);
        B = imrotate(A,angle,varargin{:});
        B = permute(B,[1 3 2]);
    case 'y'
        A = permute(A,[2 3 1]);
        B = imrotate(A,angle,varargin{:});
        B = permute(B,[2 3 1]);
    case 'z'
        B = imrotate(A,angle,varargin{:});
end
