clear all

% create table of typical coupling values and structures
min_abundance_take_into_account=1.5;%Hz Minimal natural abundance to take into account structure
min_coupling_take_into_account=1;%Hz Minimal natural abundance to take into account structure
[J_structure,typ_coupling]=create_table_of_typical_J_coupling(min_abundance_take_into_account,min_coupling_take_into_account);

% set options
opt.validate=zeros(1,3);
opt.validate(1,1)=0;% validate 1D spectrum as whole when plotting
opt.validate(1,2)=0;% validate 2D spectrum as whole when plotting
opt.validate(1,3)=0;% validate 1d multiplets... (can be 1H but also 13C)
opt.validate(1,4)=1;% validate presence of signals in 2d (can be 1H but also 13C)
opt.validate=opt.validate*0;%%% CANCEL
opt.plot_1D=1;
opt.plot_2D=opt.plot_1D;
opt.dim1=1;
opt.dim2=2;
opt.dim3=3;
%opt.plot_2D=1;
opt.only_read_data=0;
opt.draw_verbose=1*2+0*4;%0-4;this is to plot the 3d structure of the molecule in the molblok 1-4
opt.reading_sdf_file_verbose=0*2;
opt.main_prog_verbose=1;
opt.generate_eps_files=1;

%% create and manage archive %dj version
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
min_lev=10*4;

min_coupling_for_table_Jxy_HH_coupling=0.5;%min coupling for table table_Jxy used for display coupling network
atoms_no_exchange=[6 ];%for cdcl3
average_ch2_chem_shift=0;
%atoms_no_exchange=[6 7 8];%for dmso
list_location{1,1}='/Volumes/s-chior-jeannera/group_jeannerat/pupier/mnova_conversion_work';
list_location{1,2}='/Volumes/s-chior-jeannera/group_jeannerat/pupier/mnova_conversion_work/compounds/working_directory_pupier/';
list_location{1,3}='/Volumes/san256/users_for_mac_system_macPro/jeannerat/mygit/NMReDATA/examples_of_NMRrecords_and_nmredata_sdf_files/ethanol_from_DFT_GIAO_dft/';
for super_loop_over_location=1:size(list_location,2)
    folder_location=list_location{1,super_loop_over_location};
    if ~exist(folder_location,'dir')
        folder_location='./NMR_Records_folder/';
    end
    
    dataset_name='HAP_benzo(a)pyrene_assignments';
    for looop_test_zip_sdf=[ 1]%1 ziped   0 .sdf
        if looop_test_zip_sdf==1
            F = dir([  folder_location    '*.zip'] );
        else
            F = dir([  folder_location    '*.sdf'] );
        end
        for ii = 1:length(F)
            dataset_name=F(ii).name(1:end-4);
            if opt.main_prog_verbose
                if looop_test_zip_sdf==1
                    disp(['Try to read ziped dataset: ' dataset_name])
                else
                    disp(['Try to read .sdf file' dataset_name])
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%% experimental data
            if looop_test_zip_sdf==1
                folder_location_full=[folder_location dataset_name '/'];
                F2 = dir([  folder_location_full    'compound*.nmredata.sdf'] );
                
            else
                folder_location_full=[folder_location ];
                F2 = dir([  folder_location_full   '/' dataset_name  '.sdf'] );
                
            end
            for endloop_over_compound_of_record = 1:length(F2)
                part=F2(endloop_over_compound_of_record).name(1:end);
                
                %file_name=[folder_location_full 'compound1.nmredata.sdf'];
                file_name=[folder_location_full part];
                if opt.reading_sdf_file_verbose
                    disp(['> Reading file : ' file_name ' ' num2str(loop_over_compound_of_record) ' of ' num2str(length(F2))])
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [super_obj, returned_value, text_of_the_problem]=check_nmredata_sdf_file(file_name,opt);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if returned_value==0
                    warning(text_of_the_problem)
                end
                if isfield(super_obj{1}.structure,'nb_bond_between')
                    if opt.reading_sdf_file_verbose
                        disp(['> generate_table_of_nb_correl_labels  ' ])
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    draw_structure(super_obj{1},opt);
                    if opt.generate_eps_files
                        mkdir('./plots_of_figures')
                        if  exist('OCTAVE_VERSION', 'builtin') ~= 0
                            print(['./plots_of_figures/O_' dataset_name '_compound' num2str(endloop_over_compound_of_record) '_structure.eps'],'-color');%for octave
                        else
                            print(['./plots_of_figures/M_' dataset_name '_compound' num2str(endloop_over_compound_of_record) '_structure.eps'],'-depsc');%for matlab
                        end
                    end
                    [tma,tmb]=generate_table_of_nb_correl_labels(super_obj{1});
                    super_obj{1}.structure.nb_bond_between_atoms_including_implicit_H=tma;
                    super_obj{1}.structure.ref_from_label_number_to_nb_bond_table=tmb;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    error('problem');
                end
                
                if opt.reading_sdf_file_verbose
                    disp(['> End reading file : ' file_name ])
                end
                if ~opt.only_read_data
                    list_1d=[]; pointer_1d=1;
                    list_2d=[]; pointer_2d=1;
                    %% plot spectra
                    for loo=3:size(super_obj,2)%search all spectral objects (assignment is object 1 and J is object 2)
                        disp(['testing element ' num2str(loo) ' ' super_obj{loo}.tag_name ]);
                        
                        if isfield(super_obj{loo},'spectrum_location')
                            actualfile=super_obj{loo}.spectrum_location;
                            actualfile=actualfile(6:end);
                            % disp(['found field ' actualfile ]);
                            OK=0;
                            if exist([folder_location_full actualfile],'dir')
                                %disp(['found spectrum directory ' [folder_location_full actualfile] ]);
                                if isfield(super_obj{loo},'nb_dim')
                                    %  disp(['found nb_dim ' super_obj{loo}.nb_dim]);
                                    pieces=strsplit(actualfile,'/');exp_name=pieces{1,1};expno=str2num(pieces{1,2});procno=str2num(pieces{1,4});
                                    if super_obj{loo}.nb_dim == 1
                                        
                                        if exist([folder_location_full actualfile '/1r'],'file')
                                            disp(['Plotting spectrum ' num2str(loo) ' ' super_obj{loo}.tag_name ' Folder ' super_obj{loo}.spectrum_location ' exists with 1r file']);
                                            super_obj{loo}.spectrum=read_data_bruker([folder_location_full exp_name '/'],expno,procno,0);
                                            figure(loo);clf;
                                            if opt.plot_1D
                                                plot_1d(super_obj{loo}.spectrum,0,[],[],loo);
                                            end
                                            title(super_obj{loo}.tag_name,'interpreter','none')
                                            list_1d(pointer_1d,1)=loo;
                                            pointer_1d=pointer_1d+1;
                                            OK=1;
                                            %  answer = inputdlg('Is this spectrum a +D ?1h? spectrum','Check 1D spectrum',[1 40],defaultans,options);
                                            if opt.validate(1,1)>0
                                                choice = questdlg(['Is this spectrum really a ' super_obj{loo}.tag_name ' '], 'Validation of 1D spectrum', 'OK','not OK','OK');
                                                opt.validate(1,1)=opt.validate(1,1)-1;
                                            end
                                        else
                                            disp(['No 1r file found in folder ' super_obj{loo}.spectrum_location ]);
                                        end
                                    end
                                    if super_obj{loo}.nb_dim == 2
                                        if exist([folder_location_full actualfile '/2rr'],'file')
                                            disp(['Plotting spectrum ' num2str(loo) ' ' super_obj{loo}.tag_name ' Folder ' super_obj{loo}.spectrum_location ' exists with 2rr file']);
                                            super_obj{loo}.spectrum=read_data_bruker([folder_location_full exp_name '/'],expno,procno,0);
                                            [super_obj{loo}.spectrum.noise_level, super_obj{loo}.spectrum.list_peaks]=determine_noise_level(super_obj{loo}.spectrum);
                                            
                                            super_obj{loo}.spectrum.cont_level_list=generate_contour_level_list(super_obj{loo}.spectrum.noise_level * min_lev ,max(max(super_obj{loo}.spectrum.list_peaks)) , 2.00);
                                            figure(loo);clf;
                                            if opt.plot_2D
                                                plot_2d_interp(super_obj{loo}.spectrum,0,[],[],loo);
                                            end
                                            title(super_obj{loo}.tag_name,'interpreter','none')
                                            OK=1;
                                            list_2d(pointer_2d,1)=loo;
                                            pointer_2d=pointer_2d+1;
                                            if opt.validate(1,2)>0
                                                choice = questdlg(['Is this spectrum really a ' super_obj{loo}.tag_name ' '], 'Validation of 2D spectrum', 'OK','not OK','OK');
                                                opt.validate(1,2)=opt.validate(1,2)-1;
                                                
                                            end
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
                    
                    %% plot NMReDATA
                    %% plot 1D NMReDATA
                    used_color='k';
                    for loop_over_spectra=1:pointer_1d-1
                        
                        plot_1D_spectrum_and_display_data(super_obj{list_1d(loop_over_spectra,1)},list_1d(loop_over_spectra,1),used_color,opt);
                        if opt.generate_eps_files
                            mkdir('./plots_of_figures')
                            if  exist('OCTAVE_VERSION', 'builtin') ~= 0
                                print(['./plots_of_figures/O_' dataset_name '_compound' num2str(endloop_over_compound_of_record) '_'  super_obj{list_1d(loop_over_spectra,1)}.tag_name '.eps'],'-color');%for octave
                            else
                                print(['./plots_of_figures/M_' dataset_name '_compound' num2str(endloop_over_compound_of_record) '_'  super_obj{list_1d(loop_over_spectra,1)}.tag_name '.eps'],'-depsc');%for matlab
                            end
                        end
                    end
                    %% plot 2D NMReDATA
                    for loop_over_spectra=1:pointer_2d-1
                        
                        plot_2D_spectrum_and_display_data(super_obj,list_2d(loop_over_spectra,1),used_color,opt);
                        if opt.generate_eps_files
                            mkdir('./plots_of_figures')
                            if  exist('OCTAVE_VERSION', 'builtin') ~= 0
                                print(['./plots_of_figures/O_' dataset_name '_compound' num2str(endloop_over_compound_of_record) '_'  super_obj{list_2d(loop_over_spectra,1)}.tag_name '.eps'],'-color');%for octave
                            else
                                print(['./plots_of_figures/M_' dataset_name '_compound' num2str(endloop_over_compound_of_record) '_'  super_obj{list_2d(loop_over_spectra,1)}.tag_name '.eps'],'-depsc');%for matlab
                            end
                        end
                    end
                end
            end
        end
    end
end