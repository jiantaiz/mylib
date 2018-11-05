function [D, Cx, Cy, Cz, OSS] = spua_curl_slim(xdata, ydata, zdata, fovs, spua_flag, flag_2d, S_kernel,td_mode,S_slim)
%SPUA_CURL_SLIM -- Calculate the (time-domain or 1st harmonic) curl and divergence of MRE phase data + OSS (SLIM support).
%
%[D, Cx, Cy, Cz, OSS] = spua_curl_slim(X, Y, Z, fovs, spua_flag, flag_2d, S_kernel,td_mode,S_slim)
%  This function calculates the (time-domain or 1st harmonic) curl and divergence of MRE phase data.
%  This function also returns the octahedral shear strain over time.
%  The components of the MRE vector motion are input as "X", "Y", and "Z"
%  and are treated as Ny x Nx x Nz x Nt real-valued phase data.
%  "fovs" is a 3-element vector containing the FOV's in the X, Y, and Z
%  directions, in that order.
%    If "spua_flag" is 1, then the phase data are converted to complex values
%  and processed as if trying to do the calculations without unwrapping the phase.
%    If "spua_flag" is 0, then the phase data are left as real-valued data
%  and the necessary derivatives for the calculation are calculated directly
%  on the input phase data, so residual phase wraps will cause discontinuities
%  in the results.
%    If "flag_2d" is 1, then derivatives in the Z direction are not performed
%  and are assumed to be zero, plus the derivative kernels are 2D.
%    If S_kernel is absent or empty, then the default central-difference derivative
%  calculations will be performed.  Otherwise, it will use S_kernel.GRx/GRy/GRz
%  as the convolution kernels.
%    td_mode (time domain mode) = 0 (default) => output first harmonic results
%    td_mode = 1 => output time domain of full curl and divergence signal
%    td_mode = 2 => output time domain version of 1st harmonic signals
%    If the "S_slim.slim_flag" is 1, then SLIM MRE data is assumed.  Thus the X, Y, and Z
%  data are simultaneously encoded as different harmonics in xdata, and ydata and zdata
%  are empty.  The derivates for the curl have to be calculated first on this
%  composite data and then separated via the FT.
%    "S_slim.xyz" = [2,3,4], for example, to set X, Y, Z as the 2nd, 3rd, and 4th
%  harmonics of the original data.
%    "S_slim.xyz_amp" = [1,1,-1], for example, to flip the polarity of Z.
%
%  NOTE: THIS SUBROUTINE MAKES NO EFFORT TO INSURE THAT THE INPUT DATA
%  ARE IN AN APPROPRIATE RIGHT-HANDED COORDINATE SYSTEM FOR THE RESULT TO
%  TO ACTUALLY MAKE SENSE OR BE CORRECT.  THEREFORE, THE USER IS RESPONSIBLE
%  FOR CHECKING TO SEE IF THE DATA NEED TO BE SWAPPED, TRANSPOSED, AND/OR
%  NEGATED TO OBTAIN A RIGHT-HANDED COORDINATE SYSTEM.  THIS IS TYPICALLY DONE
%  BY TRYING DIFFERENCE COMBINATIONS, ORIENTATIONS, AND POLARIZATIONS OF THE
%  INPUT DATA AND LOOKING AT THE CALCULATED DIVERGENCE DATA.  THE DIVERGENCE
%  OF THE DISPLACEMENT DATA IS TYPICALLY SMALL AND UNIFORM.  IF ONE SEES EVIDENCE
%  OF WAVES THAT LOOK LIKE THE CX, CY, AND CZ DATA, THEN THIS IS AN INDICATION
%  THAT THE DATA WERE NOT IN THE CORRECT COORDINATE SYSTEM.
%
%  Kevin Glaser -- June 1, 2010
%    9/7/2011 -- Added 2D flag
%    5/25/2012 -- Added OSS
%    6/20/2013 -- Added arbitrary gradient kernel and td_mode
%

dbg_flag = 0; D = []; Cx = []; Cy = []; Cz = []; OSS = [];

% Set default values for variables not input into function.
if (nargin < 9)
    S_slim = [];
    if (nargin < 8)
        td_mode = 0;
        if (nargin < 7)
            S_kernel = [];
            if (nargin < 6)
                flag_2d = 0;
                if (nargin < 5)
                    spua_flag = 1;
                    if (nargin < 4)
                        fovs = [0.2, 0.2, 0.064];
                        if (nargin < 3)
                            zdata = [];
                            if (nargin < 2);
                                ydata = [];
                                if (nargin < 1);
                                    xdata = [];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% Figure out image size based on which displacements have been provided.
if (isempty(xdata))
    if (isempty(ydata))
        if (isempty(zdata))
            error('SPUA_CURL: All input datasets are empty.');
        else
            Nx = size(zdata,2); Ny = size(zdata,1); Nz = size(zdata,3); Nt = size(zdata,4);
        end
    else
        Nx = size(ydata,2); Ny = size(ydata,1); Nz = size(ydata,3); Nt = size(ydata,4);
    end
else
    Nx = size(xdata,2); Ny = size(xdata,1); Nz = size(xdata,3); Nt = size(xdata,4);
end

if (isempty(xdata)), xdata = repmat(single(0),[Ny,Nx,Nz,Nt]); end
if (isempty(ydata)), ydata = repmat(single(0),[Ny,Nx,Nz,Nt]); end
if (isempty(zdata)), zdata = repmat(single(0),[Ny,Nx,Nz,Nt]); end

xdata = single(xdata); ydata = single(ydata); zdata = single(zdata);

fovx = fovs(1); fovy = fovs(2); fovz = fovs(3);
dx = fovx/Nx; dy = fovy/Ny; dz = fovz/Nz;

if (spua_flag)  % data must be made complex
    if (dbg_flag==1), disp('. . Using SPUA processing'); end
    xdata = exp(i*xdata); ydata = exp(i*ydata); zdata = exp(i*zdata); uw_flag = 0;
else  % data are left real
    if (dbg_flag==1), disp('. . Using non-SPUA processing'); end
    uw_flag = 1;
end

if (isempty(S_slim))
    S_slim.slim_flag = 0; S_slim.xyz = [2,3,4]; S_slim.xyz_amp = [1,1,1];
end
slim_flag = S_slim.slim_flag; slim_xyz = S_slim.xyz; slim_xyz_amp = S_slim.xyz_amp;

% Create first-derivative convolution kernels
if (isempty(S_kernel))
    if (flag_2d == 1)
        base_gd = [-0.5, 0, 0.5];
        GRx = repmat(base_gd/(dx),[3,1])/(3);
        GRy = repmat(reshape(base_gd/(dy),[3,1]),[1,3])/(3);
        GRz = repmat(single(0),[3,3]);
    else
        base_gd = [-0.5, 0, 0.5];
        GRx = repmat(base_gd/(dx),[3,1,3])/(3^2);
        GRy = repmat(reshape(base_gd/(dy),[3,1,1]),[1,3,3])/(3^2);
        GRz = repmat(reshape(base_gd/(dz),[1,1,3]),[3,3,1])/(3^2);
    end
else
    GRx = S_kernel.GRx; GRy = S_kernel.GRy; GRz = S_kernel.GRz;
end


if (Nt > 2)
    % Calculate the divergence of the displacement data
    P=xdata; if (uw_flag); Tx_x = -convn(P,GRx,'same'); else Tx_x = real(-i*conj(P).*(-convn(P,GRx,'same'))); end, ft_phs = ifft(Tx_x,[],4); Ux_x = ft_phs(:,:,:,2);
    if (slim_flag == 1)
        Ux_x = slim_xyz_amp(1)*ft_phs(:,:,:,slim_xyz(1)); Uy_x = slim_xyz_amp(2)*ft_phs(:,:,:,slim_xyz(2)); Uz_x = slim_xyz_amp(3)*ft_phs(:,:,:,slim_xyz(3));
        if (uw_flag); Tx_y = -convn(P,GRy,'same'); else Tx_y = real(-i*conj(P).*(-convn(P,GRy,'same'))); end, ft_phs = ifft(Tx_y,[],4); Ux_y = ft_phs(:,:,:,2);
        Ux_y = slim_xyz_amp(1)*ft_phs(:,:,:,slim_xyz(1)); Uy_y = slim_xyz_amp(2)*ft_phs(:,:,:,slim_xyz(2)); Uz_y = slim_xyz_amp(3)*ft_phs(:,:,:,slim_xyz(3));
        if (uw_flag); Tx_z = -convn(P,GRz,'same'); else Tx_z = real(-i*conj(P).*(-convn(P,GRz,'same'))); end, ft_phs = ifft(Tx_z,[],4); Ux_z = ft_phs(:,:,:,2);
        Ux_z = slim_xyz_amp(1)*ft_phs(:,:,:,slim_xyz(1)); Uy_z = slim_xyz_amp(2)*ft_phs(:,:,:,slim_xyz(2)); Uz_z = slim_xyz_amp(3)*ft_phs(:,:,:,slim_xyz(3));
        D = Ux_x + Uy_y + Uz_z; if (nargout == 1), return; end
        Cx = Uz_y - Uy_z; if (nargout == 2), return; end
        Cy = Ux_z - Uz_x; if (nargout == 3), return; end
        Cz = Uy_x - Ux_y; if (nargout == 4), return; end
        OSS = []; return;
    end
    
    if (td_mode == 1), Ux_x = Tx_x; elseif (td_mode == 2), ft_phs(:,:,:,setdiff(1:Nt,[2,Nt])) = 0; Ux_x = real(fft(ft_phs,[],4)); end
    P=ydata; if (uw_flag); Ty_y = -convn(P,GRy,'same'); else Ty_y = real(-i*conj(P).*(-convn(P,GRy,'same'))); end, ft_phs = ifft(Ty_y,[],4); Uy_y = ft_phs(:,:,:,2);
    if (td_mode == 1), Uy_y = Ty_y; elseif (td_mode == 2), ft_phs(:,:,:,setdiff(1:Nt,[2,Nt])) = 0; Uy_y = real(fft(ft_phs,[],4)); end
    
    if (flag_2d == 1)
        Uz_z = repmat(single(0),[Ny,Nx,Nz]); Tz_z = repmat(single(0),size(zdata));
        if ((td_mode==1)||(td_mode==2)), Uz_z = Tz_z; end
    else
        P=zdata; if (uw_flag); Tz_z = -convn(P,GRz,'same'); else Tz_z = real(-i*conj(P).*(-convn(P,GRz,'same'))); end, ft_phs = ifft(Tz_z,[],4); Uz_z = ft_phs(:,:,:,2);
        if (td_mode == 1), Uz_z = Tz_z; elseif (td_mode == 2), ft_phs(:,:,:,setdiff(1:Nt,[2,Nt])) = 0; Uz_z = real(fft(ft_phs,[],4)); end
    end
    DT = Tx_x + Ty_y + Tz_z; D = Ux_x + Uy_y + Uz_z; if (nargout == 1), return; end
    %clear Ux_x Uy_y Uz_z
    
    % Calculate the curl of the displacement data
    P=zdata; if (uw_flag); Tz_y = -convn(P,GRy,'same'); else Tz_y = real(-i*conj(P).*(-convn(P,GRy,'same'))); end, ft_phs = ifft(Tz_y,[],4); Uz_y = ft_phs(:,:,:,2);
    if (td_mode == 1), Uz_y = Tz_y; elseif (td_mode == 2), ft_phs(:,:,:,setdiff(1:Nt,[2,Nt])) = 0; Uz_y = real(fft(ft_phs,[],4)); end
    
    if (flag_2d == 1)
        Uy_z = repmat(single(0),[Ny,Nx,Nz]); Ty_z = repmat(single(0),size(ydata));
        if ((td_mode==1)||(td_mode==2)), Uy_z = Ty_z; end
    else
        P=ydata; if (uw_flag); Ty_z = -convn(P,GRz,'same'); else Ty_z = real(-i*conj(P).*(-convn(P,GRz,'same'))); end, ft_phs = ifft(Ty_z,[],4); Uy_z = ft_phs(:,:,:,2);
        if (td_mode == 1), Uy_z = Ty_z; elseif (td_mode == 2), ft_phs(:,:,:,setdiff(1:Nt,[2,Nt])) = 0; Uy_z = real(fft(ft_phs,[],4)); end
    end
    CTx = Tz_y - Ty_z; Cx = Uz_y - Uy_z; if (nargout == 2), return; end
    %clear Uz_y Uy_z
    
    if (flag_2d == 1)
        Ux_z = repmat(single(0),[Ny,Nx,Nz]); Tx_z = repmat(single(0),size(xdata));
        if ((td_mode==1)||(td_mode==2)), Ux_z = Tx_z; end
    else
        P=xdata; if (uw_flag); Tx_z = -convn(P,GRz,'same'); else Tx_z = real(-i*conj(P).*(-convn(P,GRz,'same'))); end, ft_phs = ifft(Tx_z,[],4); Ux_z = ft_phs(:,:,:,2);
        if (td_mode == 1), Ux_z = Tx_z; elseif (td_mode == 2), ft_phs(:,:,:,setdiff(1:Nt,[2,Nt])) = 0; Ux_z = real(fft(ft_phs,[],4)); end
    end
    P=zdata; if (uw_flag); Tz_x = -convn(P,GRx,'same'); else Tz_x = real(-i*conj(P).*(-convn(P,GRx,'same'))); end, ft_phs = ifft(Tz_x,[],4); Uz_x = ft_phs(:,:,:,2);
    if (td_mode == 1), Uz_x = Tz_x; elseif (td_mode == 2), ft_phs(:,:,:,setdiff(1:Nt,[2,Nt])) = 0; Uz_x = real(fft(ft_phs,[],4)); end
    CTy = Tx_z - Tz_x; Cy = Ux_z - Uz_x; if (nargout == 3), return; end
    %clear Ux_z Uz_x
    
    P=xdata; if (uw_flag); Tx_y = -convn(P,GRy,'same'); else Tx_y = real(-i*conj(P).*(-convn(P,GRy,'same'))); end, ft_phs = ifft(Tx_y,[],4); Ux_y = ft_phs(:,:,:,2);
    if (td_mode == 1), Ux_y = Tx_y; elseif (td_mode == 2), ft_phs(:,:,:,setdiff(1:Nt,[2,Nt])) = 0; Ux_y = real(fft(ft_phs,[],4)); end
    
    P=ydata; if (uw_flag); Ty_x = -convn(P,GRx,'same'); else Ty_x = real(-i*conj(P).*(-convn(P,GRx,'same'))); end, ft_phs = ifft(Ty_x,[],4); Uy_x = ft_phs(:,:,:,2);
    if (td_mode == 1), Uy_x = Ty_x; elseif (td_mode == 2), ft_phs(:,:,:,setdiff(1:Nt,[2,Nt])) = 0; Uy_x = real(fft(ft_phs,[],4)); end
    
    CTz = Ty_x - Tx_y; Cz = Uy_x - Ux_y; if (nargout == 4), return; end
    %clear Ux_y Uy_x
    OSS = (2/3)*sqrt((Tx_x-Ty_y).^2 + (Tx_x-Tz_z).^2 + (Ty_y-Tz_z).^2 + 6*((0.5*(Tx_y+Ty_x)).^2 + (0.5*(Tx_z+Tz_x)).^2 + (0.5*(Ty_z+Tz_y)).^2));
elseif (Nt == 2)
    % Calculate the divergence of the displacement data
    P=xdata; if (uw_flag); ft_phs = -convn(P,GRx,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRx,'same'))); end, Ux_x = ft_phs(:,:,:,1)+i*ft_phs(:,:,:,2); Tx_x = ft_phs;
    P=ydata; if (uw_flag); ft_phs = -convn(P,GRy,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRy,'same'))); end, Uy_y = ft_phs(:,:,:,1)+i*ft_phs(:,:,:,2); Ty_y = ft_phs;
    if (flag_2d == 1)
        Uz_z = repmat(single(0),[Ny,Nx,Nz]); Tz_z = repmat(single(0),size(zdata));
    else
        P=zdata; if (uw_flag); ft_phs = -convn(P,GRz,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRz,'same'))); end, Uz_z = ft_phs(:,:,:,1)+i*ft_phs(:,:,:,2); Tz_z = ft_phs;
    end
    DT = Tx_x + Ty_y + Tz_z; D = Ux_x + Uy_y + Uz_z; if (nargout == 1), return; end
    %clear Ux_x Uy_y Uz_z

    % Calculate the curl of the displacement data
    P=zdata; if (uw_flag); ft_phs = -convn(P,GRy,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRy,'same'))); end, Uz_y = ft_phs(:,:,:,1)+i*ft_phs(:,:,:,2); Tz_y = ft_phs;
    if (flag_2d == 1)
        Uy_z = repmat(single(0),[Ny,Nx,Nz]); Ty_z = repmat(single(0),size(ydata));
    else
        P=ydata; if (uw_flag); ft_phs = -convn(P,GRz,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRz,'same'))); end, Uy_z = ft_phs(:,:,:,1)+i*ft_phs(:,:,:,2); Ty_z = ft_phs;
    end
    CTx = Tz_y - Ty_z; Cx = Uz_y - Uy_z; if (nargout == 2), return; end
    %clear Uz_y Uy_z
    if (flag_2d == 1)
        Ux_z = repmat(single(0),[Ny,Nx,Nz]); Tx_z = repmat(single(0),size(xdata));
    else
        P=xdata; if (uw_flag); ft_phs = -convn(P,GRz,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRz,'same'))); end, Ux_z = ft_phs(:,:,:,1)+i*ft_phs(:,:,:,2); Tx_z = ft_phs;
    end
    P=zdata; if (uw_flag); ft_phs = -convn(P,GRx,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRx,'same'))); end, Uz_x = ft_phs(:,:,:,1)+i*ft_phs(:,:,:,2); Tz_x = ft_phs;
    CTy = Tx_z - Tz_x; Cy = Ux_z - Uz_x; if (nargout == 3), return; end    
    %clear Ux_z Uz_x
    P=xdata; if (uw_flag); ft_phs = -convn(P,GRy,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRy,'same'))); end, Ux_y = ft_phs(:,:,:,1)+i*ft_phs(:,:,:,2); Tx_y = ft_phs;
    P=ydata; if (uw_flag); ft_phs = -convn(P,GRx,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRx,'same'))); end, Uy_x = ft_phs(:,:,:,1)+i*ft_phs(:,:,:,2); Ty_x = ft_phs;
    CTz = Ty_x - Tx_y; Cz = Uy_x-Ux_y; if (nargout == 4), return; end
    OSS = (2/3)*sqrt((Tx_x-Ty_y).^2 + (Tx_x-Tz_z).^2 + (Ty_y-Tz_z).^2 + 6*((0.5*(Tx_y+Ty_x)).^2 + (0.5*(Tx_z+Tz_x)).^2 + (0.5*(Ty_z+Tz_y)).^2));
elseif (Nt == 1)
    % Calculate the divergence of the displacement data
    P=xdata; if (uw_flag); ft_phs = -convn(P,GRx,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRx,'same'))); end, Ux_x = ft_phs;
    P=ydata; if (uw_flag); ft_phs = -convn(P,GRy,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRy,'same'))); end, Uy_y = ft_phs;
    if (flag_2d == 1)
        Uz_z = repmat(single(0),[Ny,Nx,Nz]);
    else
        P=zdata; if (uw_flag); ft_phs = -convn(P,GRz,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRz,'same'))); end, Uz_z = ft_phs;
    end
    Tx_x = Ux_x; Ty_y = Uy_y; Tz_z = Uz_z; D = Ux_x + Uy_y + Uz_z; DT = D; if (nargout == 1), return; end
    %clear Ux_x Uy_y Uz_z

    % Calculate the curl of the displacement data
    P=zdata; if (uw_flag); ft_phs = -convn(P,GRy,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRy,'same'))); end, Uz_y = ft_phs;
    if (flag_2d == 1)
        Uy_z = repmat(single(0),[Ny,Nx,Nz]);
    else
        P=ydata; if (uw_flag); ft_phs = -convn(P,GRz,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRz,'same'))); end, Uy_z = ft_phs;
    end
    Tz_y = Uz_y; Ty_z = Uy_z; Cx = Uz_y-Uy_z; CTx = Cx; if (nargout == 2), return; end
    %clear Uz_y Uy_z
    if (flag_2d == 1)
        Ux_z = repmat(single(0),[Ny,Nx,Nz]);
    else
        P=xdata; if (uw_flag); ft_phs = -convn(P,GRz,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRz,'same'))); end, Ux_z = ft_phs;
    end
    P=zdata; if (uw_flag); ft_phs = -convn(P,GRx,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRx,'same'))); end, Uz_x = ft_phs;
    Tx_z = Ux_z; Tz_x = Uz_x; Cy = Ux_z-Uz_x; CTy = Cy; if (nargout == 3), return; end    
    %clear Ux_z Uz_x
    P=xdata; if (uw_flag); ft_phs = -convn(P,GRy,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRy,'same'))); end, Ux_y = ft_phs;
    P=ydata; if (uw_flag); ft_phs = -convn(P,GRx,'same'); else ft_phs = real(-i*conj(P).*(-convn(P,GRx,'same'))); end, Uy_x = ft_phs;
    Tx_y = Ux_y; Ty_x = Uy_x; Cz = Uy_x-Ux_y; CTz = Cz; return;
    OSS = (2/3)*sqrt((Tx_x-Ty_y).^2 + (Tx_x-Tz_z).^2 + (Ty_y-Tz_z).^2 + 6*((0.5*(Tx_y+Ty_x)).^2 + (0.5*(Tx_z+Tz_x)).^2 + (0.5*(Ty_z+Tz_y)).^2));
end



