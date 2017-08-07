function [system]=read_dft_nmr_data_fn(dft_folder)
%% this function reads data from files in which gaussian data were 
%% extracted from DFT/GIAO gaussian output files
%% using "report.c" 
%% compile using : gcc -o report report.c
%% call from the folder where the .log files is/are
%% for example :  ../report *nmrJ.log
%% it will (among others) the free files read here:
%% log.txtj 
%% log.txt 
%% log.min.CXYt 


file=[ dft_folder 'log.txtj'];
if ~exist(file,'file')
    disp(['file ' file ' Not found !'])
end
nmrJ_log=load(file);
file=[ dft_folder 'log.txt'];
if ~exist(file,'file')
    disp(['file ' file ' Not found !'])
end
fid=fopen(file);
C = textscan(fid, '%d%s%f');
fclose(fid);
%num=C{1};
%atom=C{2};
cs=C{3};%chemical shifts

file=[ dft_folder 'log.min.XYZt'];
if ~exist(file,'file')
    disp(['file ' file ' Not found !'])
end
fid=fopen(file);
C = textscan(fid, '%d%f%f%f');
fclose(fid);
xyz=[ C{2} C{3} C{4}];
atom_number=C{1} ;
system.xyz=xyz;
system.atom_number=atom_number;
system.J=nmrJ_log;
system.cs=cs;

end