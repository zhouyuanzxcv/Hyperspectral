% This script will compile all the C files of the registration methods
cd('functions_nonrigid');
files=dir('*.c');
for i=1:length(files)
    filename=[files(i).name];
    if(length(filename)<19||(~strcmpi(filename(1:19),'image_interpolation')))
        disp(['compiling : ' filename]);
        mex(filename,'image_interpolation.c','-v');
    end
end
cd('..');

cd('functions_affine');
files=dir('*.c');
for i=1:length(files)
    filename=[files(i).name];
    if(length(filename)<19||(~strcmpi(filename(1:19),'image_interpolation')))
        disp(['compiling : ' filename]);
        mex(filename,'image_interpolation.c','-v');
    end
end
cd('..');

cd('functions')
files=dir('*.c');
for i=1:length(files)
    filename=[files(i).name];
    disp(['compiling : ' filename]);
    mex(filename,'-v');
end
cd('..');


