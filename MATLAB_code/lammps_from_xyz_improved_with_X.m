function lammps_from_xyz_improved_with_X(fname, natoms);
% Convert standard XYZ file to LAMMPS input file format
% Example - 
%                  input file is:                  fname='min_end'; (Omit the ending .xyz when calling this function - I add it inside this code)
%                  ouput file is:                 fname='min_end_lmps.in';
%                  calling this function:     lammps_from_xyz_improved_with_X('min_end', 100)   *This is an example with 100 atoms in the xyz file

dir_out='';
dir_in='';
fname_in=[dir_in,fname,'.xyz'];
fname_out=[dir_out,fname,'_lmps.in'];
[x,y,z,atom_name]=read_data(fname_in);
check_size=natoms-length(x);
if(check_size~=0)
    disp(['problem with number of atoms in input - should be = ', num2str(length(x))])
    pause
end
xmin=min(x);
xmax=max(x);
ymin=min(y);
ymax=max(y);
zmin=min(z);
zmax=max(z);
disp(['xmin-xmax = ',num2str(xmin),' - ', num2str(xmax)])
disp(['ymin-ymax = ',num2str(ymin),' - ', num2str(ymax)])
disp(['zmin-zmax = ',num2str(zmin),' - ', num2str(zmax)])

atom_id=linspace(1,natoms,natoms);
atom_0=zeros(natoms,1);
for k=1:length(x)
    switch atom_name(k)
        case 'C'
            atom_num(k)=1;
        case 'H'
            atom_num(k)=2;
        case 'O'
            atom_num(k)=3;
        case 'N'
            atom_num(k)=4;
        case 'S'
            atom_num(k)=5;
    end
end
if(length(x)==0)
    disp('forgot to erase the 2nd line in the .xyz file - ERASE it (leave empty) and try again')
end

fid = fopen(fname_out, 'wt');
for i=1:natoms
    fprintf(fid, '%d %d %d %12.6f %12.6f %12.6f \n',atom_id(i),atom_num(i),atom_0(i),x(i),y(i),z(i));
end
fclose (fid);
return

