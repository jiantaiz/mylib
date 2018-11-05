function data = load_nd_data(file,yxz,dtype,bs,offset,c)
%LOAD_ND_DATA -- Read in a N-D binary file into an array
%
%data = load_nd_data(file,yxz,dtype,bs,offset,c);
%  Read in a N-D binary file whose name is given in the string FILE.
%  The program will skip the first OFFSET number of bytes, and then
%  read in the next elements as if they were of type DTYPE (string).
%  YXZ is a n-element array giving the y, x, z, ... dimensions, respectively,
%  of the array into which the data will be read.  The resulting array
%  will be YXZ(1) rows by YXZ(2) columns by YXZ(3) "slices" ... .  BS is
%  a string indicating what kind of byte-order to use.  If BS is 'b',
%  then the 'big-endian' (standard UNIX) format will be used.  If
%  BS is 'l', then the 'little-endian' (Windows PC) format will be used.
%  C is a flag identifying the data as complex or not ('1' for complex,
%  '0' for real).  If the data are complex, it will assume that the real
%  and imaginary parts are alternating (real, imag, real, imag, . .).
%
%Defaults: OFFSET = 0, BS = 'b', DTYPE = 'uchar', YXZ = [256,256,1], C = 0.
%Note: This program assumes that the binary file is organized (e.g. in 3D) as follows
%  slice 1: row 1, row 2, . . . slice 2: row 1, row 2, . . .
%  and that Matlab reads in data as a COLUMN each time.  Therefore,
%  the program will read in row 1 of the data initially as column 1 of
%  the results array, and then later transpose the data to make column 1
%  become row 1 in the output.
%
%Kevin Glaser -- June 7, 2006
%

if (nargin < 6)
   c = 0;
   if (nargin < 5)
      offset = 0;
      if (nargin < 4)
         bs = 'b';
         if (nargin < 3)
            dtype = 'uchar';
            if (nargin < 2)
               yxz = [256,256,1];
               if (nargin < 1)
                  error('READ_BINARY_FILE: At least 1 argument required');
               end
            end
         end
      end
   end
end

fid = fopen(file,'r',bs);
if (fid == -1)
   disp(['READ_BINARY_FILE: File (',file,') was not opened correctly']);
end

status = fseek(fid,offset,-1);  % skip OFFSET bytes into file
if (status == -1)
   error(['READ_BINARY_FILE: Fseek error (offset = ',num2str(offset),')']);
end

if (length(yxz) < 3)
   yxz = [yxz,ones(1,3-length(yxz))];
end

if (prod(yxz)==0)
   error('READ_BINARY_FILE: Array dimensions cannot have 0 size');
end

xyz = yxz; xyz(1)=yxz(2); xyz(2)=yxz(1);
p_a = 1:length(xyz); p_a(1)=2; p_a(2)=1;
n_elements = prod(yxz);  % total number of elements to read in
if (c ~= 1)  % if the data are real
   [tmp, count] = fread(fid,n_elements,dtype);
   fclose('all');
   if (count ~= n_elements)
      error(['READ_BINARY_FILE: Error reading in (real) data (count,n_elements) = (',...
          num2str(count),',',num2str(n_elements),')']);
   end
   data = permute(reshape(tmp,xyz),p_a);
else  % data are complex
   [tmp, count] = fread(fid,2*n_elements,dtype);
   fclose('all');
   if (count ~= 2*n_elements)
      error(['READ_BINARY_FILE: Error reading in (complex) data (count,2*n_elements) = (',...
          num2str(count),',',num2str(2*n_elements),')']);
   end
   data = permute(reshape(tmp(1:2:end)+i*tmp(2:2:end),xyz),p_a);
end

