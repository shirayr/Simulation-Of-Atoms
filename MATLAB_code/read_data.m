function [x,y,z,atom_name]=read_data(fname)
fid = fopen(fname);
num_atoms = textscan(fid, '%f64', 1);
data=textscan(fid,'%c %f64 %f64 %f64 ');
x=data{:,2};
y=data{:,3};
z=data{:,4};
atom_name=data{:,1};
fclose(fid);