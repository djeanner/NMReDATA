clear all

%% create and manage archive
create_archive=0;% should always do this to keep the matlab files with the data generated...

if create_archive
    pack_vers='1.0';package_name=[mfilename '_package_folder_' pack_vers '/'];
    location_package=['./' package_name];
    if ~exist(location_package,'dir')
        disp('copy the following functions in the archive folder:')
        mkdir (location_package);
        [fList,pList] = matlab.codetools.requiredFilesAndProducts(mfilename);
        for i=1:size(fList,2)
            disp(fList{1,i});
            copyfile(which(fList{1,i}),location_package)
        end
    else
        [list_out]=dir(['./' mfilename '_package_folder_' pack_vers '/*.m']);
        fList=list_out;
        disp('Update the following functions in the archive folder:')
        for i=1:size(fList)
            disp(fList(i).name);
            copyfile(which(fList(i).name),location_package)
        end
    end
end

%
min_lev=10;

%param of contour plot on 1D spectru,
cont.plot_contour=2;%1 for vertical multiplet 2 for spectrum
cont.vbas=0.8;
cont.vwidth=0.2;

%param for J coupling representations
moulin_dft.from_y=1.05;%position on paper
moulin_dft.width_y=-0.2;%if positive goes down...

% params for representation of coupling fro dft if available
moulin_dft.from_y=0;%position on paper
moulin_dft.width_y=-23;%if positive goes down...
moulin_dft.ticks_at=5;%position of ticks every ... Hz
moulin_dft.max_cou=23;% max coupling considered (for scales of )
moulin_dft.my_color_map=colormap(hsv(128));

dft_folder='./';
top_page=max([(moulin_dft.from_y-moulin_dft.width_y) (moulin_dft.from_y-moulin_dft.width_y) (cont.vbas+cont.vwidth)])+0.1;

min_coupling_for_table_Jxy_HH_coupling=0.5;%min coupling for table table_Jxy used for display coupling network
atoms_no_exchange=[6 ];%for cdcl3
average_ch2_chem_shift=0;
%atoms_no_exchange=[6 7 8];%for dmso

%%%%%%%%%%%%%%%%%%%%%%%% experimental data
folder_location='./NMR_Records_folder/';
folder_location_full=[folder_location 'HAP_benzo(a)pyrene_assignments/'];
file_name=[folder_location_full 'compound1.nmredata.sdf'];

[super_obj, returned_value, text_of_the_problem]=check_nmredata_sdf_file(file_name);

%% plot spectra
for loo=3:size(super_obj,2)%search all spectral objects (assignment is object 1 and J is object 2)
    if isfield(super_obj{loo},'spectrum_location')
        %  disp( super_obj{loo}.spectrum_location);
        actualfile=super_obj{loo}.spectrum_location;
        actualfile=actualfile(6:end);
        OK=0;
        if exist([folder_location_full actualfile],'dir')
            if isfield(super_obj{loo},'nb_dim')
                pieces=strsplit(actualfile,'/');exp_name=pieces{1,1};expno=str2num(pieces{1,2});procno=str2num(pieces{1,4});
                
                if super_obj{loo}.nb_dim == 1
                    
                    if exist([folder_location_full actualfile '/1r'],'file')
                        disp(['Element ' num2str(loo) ' Folder ' super_obj{loo}.spectrum_location ' exists with 1r file']);
                        super_obj{loo}.spectrum=read_data_bruker([folder_location_full exp_name '/'],expno,procno,0);
                        
                        plot_1d(super_obj{loo}.spectrum);
                        print(['./spectrum_' num2str(loo) '.eps']);%here
                        
                        OK=1;
                    else
                        disp(['No 1r file found in folder ' super_obj{loo}.spectrum_location ]);
                    end
                end
                if super_obj{loo}.nb_dim == 2
                    if exist([folder_location_full actualfile '/2rr'],'file')
                        disp(['Element ' num2str(loo) ' Folder ' super_obj{loo}.spectrum_location ' exists with 2rr file']);
                        super_obj{loo}.spectrum=read_data_bruker([folder_location_full exp_name '/'],expno,procno,0);
                        [super_obj{loo}.spectrum.noise_level, super_obj{loo}.spectrum.list_peaks]=determine_noise_level(super_obj{loo}.spectrum);
                        
                        super_obj{loo}.spectrum.cont_level_list=generate_contour_level_list(super_obj{loo}.spectrum.noise_level * min_lev ,max(max(super_obj{loo}.spectrum.list_peaks)) , 2.00);
                        plot_2d_interp(super_obj{loo}.spectrum);
                        print(['./spectrum_' num2str(loo) '.eps']);%here
                        OK=1;
                    else
                        disp(['No 2rr file found in folder ' super_obj{loo}.spectrum_location ]);
                    end
                end
            end
        else
            disp(['Element ' num2str(loo) ' Folder ' super_obj{loo}.spectrum_location ' does not exists']);
        end
    end
end




drawnow
%%%%%%%%%%%%%%%%%%%%%%%% dft data


stopsfadsfa

name='androsten';
average_ch2_chem_shift=0;
full_path='./dft_androsten/androsten_better_nmrJ.';

[system,red_system]=get_from_dft_data(full_path,min_coupling_for_table_Jxy_HH_coupling,atoms_no_exchange,average_ch2_chem_shift,name);
draw_coupling_plot_new(moulin_dft,red_system);
%draw_coupling_plot_new(moulin_dft,system);
[super_obj, returned_value, text_of_the_problem]=check_nmredata_sdf_file(['./' name '.sdf']);
super_obj{:}
[super_obj, returned_value, text_of_the_problem]=check_nmredata_sdf_file(['./' name '_implicit.sdf']);
super_obj{:}
drawnow

%ethanol
name='etoh';
average_ch2_chem_shift=1;exchange_OH_NH=1;
%  path_many_hap='/Volumes/lacie_case/work_remove_duplicate_second/folder_1_234232341234434125467856/djeanner1/xavier_cluster/gaussian/gaussian/tst/';
%  full_path=[path_many_hap 'etoh.better.nmrJ.'];%copyfile([full_path '*'],'./gaussian_files');
%  copyfile([full_path '*'],'./gaussian_files');

path_many_hap='./gaussian_files/';
full_path=[path_many_hap 'etoh.better.nmrJ.'];


[system,red_system]=get_from_dft_data(full_path,min_coupling_for_table_Jxy_HH_coupling,atoms_no_exchange,average_ch2_chem_shift,name);
draw_coupling_plot_new(moulin_dft,red_system);
check_nmredata_sdf_file(['./' name '.sdf']);
drawnow

name='cla_fig8.26';
average_ch2_chem_shift=0;
full_path=[path_many_hap 'three_better_no_oh_bondCO_bis.nmrJ.'];%copyfile([full_path '*'],'./gaussian_files');

[system,red_system]=get_from_dft_data(full_path,min_coupling_for_table_Jxy_HH_coupling,atoms_no_exchange,average_ch2_chem_shift,name);
draw_coupling_plot_new(moulin_dft,red_system);
check_nmredata_sdf_file(['./' name '.sdf']);
drawnow


name='7-12dimethylbenzaanthracene';
average_ch2_chem_shift=1;
full_path=[path_many_hap '7-12dimethylbenzaanthracene.nmrJ.'];%copyfile([full_path '*'],'./gaussian_files');

[system,red_system]=get_from_dft_data(full_path,min_coupling_for_table_Jxy_HH_coupling,atoms_no_exchange,average_ch2_chem_shift,name);
draw_coupling_plot_new(moulin_dft,red_system);
check_nmredata_sdf_file(['./' name '.sdf']);
drawnow

name='benzocchrysene';
average_ch2_chem_shift=0;
full_path=[path_many_hap 'benzocchrysene_non_flat_nmrJ.'];%copyfile([full_path '*'],'./gaussian_files');

[system,red_system]=get_from_dft_data(full_path,min_coupling_for_table_Jxy_HH_coupling,atoms_no_exchange,average_ch2_chem_shift,name);
draw_coupling_plot_new(moulin_dft,red_system);
save_sdf(system,red_system,[system.name '.sdf']);
check_nmredata_sdf_file(['./' name '.sdf']);
drawnow

name='benapyrene_nmrJ';
full_path=[path_many_hap 'benapyrene_nmrJ.'];%copyfile([full_path '*'],'./gaussian_files');

[system,red_system]=get_from_dft_data(full_path,min_coupling_for_table_Jxy_HH_coupling,atoms_no_exchange,average_ch2_chem_shift,name);
draw_coupling_plot_new(moulin_dft,red_system);
check_nmredata_sdf_file(['./' name '.sdf']);
drawnow


% process all the spectra in the gaussian_log_file_folder folder
% %
% % in_dir='./gaussian_log_file_folder/';
% % dir_out=dir(in_dir);
% %
% % for loop=1:size(dir_out,1)
% %     current_file_name=dir_out(loop).name;
% %     if contains(current_file_name,'log')
% %         if ~contains(current_file_name,'log.')
% %
% %             name=current_file_name(1:length(current_file_name)-4);
% %
% %             full_name=[in_dir  name '.'];
% %             check_needed_files_are_there=1;
% %             check_needed_files_are_there=check_needed_files_are_there*exist([full_name 'log.txtj']);
% %             if check_needed_files_are_there
% %                 disp(name)
% %
% %                 [system,red_system]=get_from_dft_data([full_name ],min_coupling_for_table_Jxy_HH_coupling,atoms_no_exchange,average_ch2_chem_shift,name);
% %                 draw_coupling_plot_new(moulin_dft,red_system);
% %             else
% %                 disp(['Not all files available for : ' [full_name 'log.txtj']])
% %
% %             end
% %             drawnow
% %             check_nmredata_sdf_file(['./' name '.sdf']);
% %             check_nmredata_sdf_file(['./' name '_implicit.sdf']);
% %         end
% %     end
% % end
