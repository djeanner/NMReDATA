clear all


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
folder_location='./nmr_spectra/benzoapyrene/';
file_name=[folder_location 'benzoapyrene.sdf'];
[super_obj, returned_value, text_of_the_problem]=check_nmredata_sdf_file(file_name);
for loo=1:size(super_obj,2)
    if isfield(super_obj{loo},'spectrum_location')
        disp( super_obj{loo}.spectrum_location);
        actualfile=super_obj{loo}.spectrum_location;
        actualfile=actualfile(6:end);
        if exist([folder_location actualfile],'dir') disp('OK folder exists');
        else disp([' folder : ' folder_location actualfile ' does not exists']) ;
        end
    end
end
drawnow
%%%%%%%%%%%%%%%%%%%%%%%% dft data

name='androsten';
average_ch2_chem_shift=0;
full_path='./dft_androsten/androsten_better_nmrJ.';
copyfile([full_path '*'],'./gaussian_files');

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
