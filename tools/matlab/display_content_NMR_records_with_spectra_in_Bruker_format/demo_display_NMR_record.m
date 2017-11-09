clear all

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
min_lev=10;

min_coupling_for_table_Jxy_HH_coupling=0.5;%min coupling for table table_Jxy used for display coupling network
atoms_no_exchange=[6 ];%for cdcl3
average_ch2_chem_shift=0;
%atoms_no_exchange=[6 7 8];%for dmso

folder_location='/Volumes/s-chior-jeannera/group_jeannerat/pupier/mnova_conversion_work/compounds/working_directory_pupier_nov.8.2017/';
folder_location='./NMR_Records_folder/';

dataset_name='HAP_benzo(a)pyrene_assignments';
F = dir([  folder_location    '*.zip'] );
for ii = 1:length(F)
   dataset_name=F(ii).name(1:end-4);
disp(['Try to read ' dataset_name])
%%%%%%%%%%%%%%%%%%%%%%%% experimental data
folder_location_full=[folder_location dataset_name '/'];
file_name=[folder_location_full 'compound1.nmredata.sdf'];

[super_obj, returned_value, text_of_the_problem]=check_nmredata_sdf_file(file_name);
list_1d=[]; pointer_1d=1;
list_2d=[]; pointer_2d=1;
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
                        figure(loo);clf;
                        
                        plot_1d(super_obj{loo}.spectrum,0,[],[],loo);
                        title(super_obj{loo}.tag_name,'interpreter','none')
                        list_1d(pointer_1d,1)=loo;
                        pointer_1d=pointer_1d+1;
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
                        figure(loo);clf;
                        
                            plot_2d_interp(super_obj{loo}.spectrum,0,[],[],loo);
                        title(super_obj{loo}.tag_name,'interpreter','none')
                        OK=1;
                        list_2d(pointer_2d,1)=loo;
                        pointer_2d=pointer_2d+1;
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
used_color='k';
for loop_over_spectra=1:pointer_1d-1
    figure(list_1d(loop_over_spectra,1));hold on;
    tmp_obj=super_obj{list_1d(loop_over_spectra,1)};
    min_chem_shift_spectrum=10000000;
    max_chem_shift_spectrum=-10000000000;
    for loop_over_peaks=1:size(tmp_obj.label,2)
        % disp(['For Peak : ' num2str(loop_over_peaks) ' Label= ' tmp_obj.label{1,loop_over_peaks}])
        %  plot(tmp_obj.chemical_shift{1,loop_over_peaks},0,[used_color '+'])
        % text(tmp_obj.chemical_shift{1,loop_over_peaks},0,tmp_obj.label{1,loop_over_peaks},'color',used_color)
        txt=tmp_obj.label{1,loop_over_peaks};
        txt=[txt ' '];
        if isfield(tmp_obj,'multiplicity')
            if tmp_obj.integral{1,loop_over_peaks}~=0
                %  txt=[txt num2str(tmp_obj.multiplicity{1,loop_over_peaks}) ' '];
            end
        end
        if isfield(tmp_obj,'integral')
            if tmp_obj.integral{1,loop_over_peaks}~=0
                %  txt=[txt 'E=' num2str(tmp_obj.integral{1,loop_over_peaks}) ' '];
            end
        end
        pos_lab=max(max(tmp_obj.spectrum.spectrum))/2;
        nb_pt_for_broadening=1+round(0.025/(max(max(tmp_obj.spectrum.scale2))-min(min(tmp_obj.spectrum.scale2)))*tmp_obj.spectrum.si2);
        %%   broadened_spectrum=conv(,ones(nb_pt_for_broadening,1),'same')/nb_pt_for_broadening;
        [del, index]=min(abs(tmp_obj.spectrum.scale2-   tmp_obj.chemical_shift{1,loop_over_peaks}));
        from=index-nb_pt_for_broadening;
        to=index+nb_pt_for_broadening;
        if from<1, from=1;end
        if to>tmp_obj.spectrum.si2, to=si2;end
        pos_lab=max(max(tmp_obj.spectrum.spectrum(from:to,1)));
        
        if isfield(tmp_obj,'intensity')
            if tmp_obj.intensity{1,loop_over_peaks}~=0
                pos_lab=tmp_obj.intensity{1,loop_over_peaks};
            end
        end
        %  disp(txt)
        list_pos=tmp_obj.chemical_shift{1,loop_over_peaks};
        
        text(list_pos,pos_lab,txt,'color',used_color)
        %% plot j structure
        if isfield(tmp_obj,'J')
            if size(tmp_obj.J,1)>=loop_over_peaks
                curJl=zeros(1,size(tmp_obj.J,2));
                for loo_j=1: size(tmp_obj.J,2)
                    totm=tmp_obj.J{loop_over_peaks,loo_j};
                    if size(totm,2)==0
                        curJl(1,loo_j)=0;
                    else
                        curJl(1,loo_j)=totm;
                    end
                end
                curJl=sort(curJl,'descend');
                for loo_j=1: size(curJl,2)
                    curJ=curJl(1,loo_j);
                    if curJ~=0
                        shiftJ=0.5*curJ/tmp_obj.spectrum.sf2;
                        alist_pos=[list_pos list_pos];
                        list_pos=[list_pos+shiftJ list_pos-shiftJ];
                        n_pos_lab=pos_lab-curJ*pos_lab*0.005;
                        for lol=1:size(alist_pos,2)
                            plot([alist_pos(1,lol) list_pos(1,lol)],[ pos_lab n_pos_lab],[used_color '-']);
                        end
                        pos_lab=n_pos_lab;
                    end
                end
                
            end
        end
        %% just little vertical line
        alist_pos=[list_pos];
                        n_pos_lab=pos_lab-4*pos_lab*0.005;
        for lol=1:size(alist_pos,2)
            plot([alist_pos(1,lol) list_pos(1,lol)],[ pos_lab n_pos_lab],[used_color '-']);
        end
        
        
        
        min_chem_shift_spectrum=min(min([min_chem_shift_spectrum list_pos]));
        max_chem_shift_spectrum=max(max([max_chem_shift_spectrum list_pos]));
    end
    margin=0.2;
    min_chem_shift_spectrum=min_chem_shift_spectrum -margin;min_chem_shift_spectrum=round(min_chem_shift_spectrum*10)/10;
    max_chem_shift_spectrum=max_chem_shift_spectrum +margin;max_chem_shift_spectrum=round(max_chem_shift_spectrum*10)/10;
    xlim([min_chem_shift_spectrum max_chem_shift_spectrum])
    
    fig_num=list_1d(loop_over_spectra,1);
     if  exist('OCTAVE_VERSION', 'builtin') ~= 0
        print(['./O_spectrum_' dataset_name '_' tmp_obj.tag_name '.eps'],'-color');%for octave
    else
        print(['./M_spectrum_' dataset_name '_' tmp_obj.tag_name '.eps'],'-depsc');%for matlab
    end
end
for loop_over_spectra=1:pointer_2d-1
    figure(list_2d(loop_over_spectra,1));hold on;
        tmp_obj=super_obj{list_2d(loop_over_spectra,1)};

    
        %%%%%%%%%%
        
        min_chem_shift_spectrum1=10000000;
        max_chem_shift_spectrum1=-10000000000;
        min_chem_shift_spectrum2=10000000;
        max_chem_shift_spectrum2=-10000000000;
        for loop_over_peaks=1:size(tmp_obj.signalf1,2)
                        list_pos1=0;
                        
                        list_pos2=0;
                        
                        % disp(['For Peak : ' num2str(loop_over_peaks) ' Label= ' tmp_obj.label{1,loop_over_peaks}])
                        %  plot(tmp_obj.chemical_shift{1,loop_over_peaks},0,[used_color '+'])
                        % text(tmp_obj.chemical_shift{1,loop_over_peaks},0,tmp_obj.label{1,loop_over_peaks},'color',used_color)
                        txt='';
                        txt=[txt tmp_obj.signalf1{1,loop_over_peaks} '/'];
                        txt=[txt tmp_obj.signalf2{1,loop_over_peaks}];
                        for loop_over_assigned_signals=1:size(super_obj{1,1}.label_signal,2)
                            if size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},2)>0
                            if size(tmp_obj.signalf1{1,loop_over_peaks},2)>0
                              %  disp(['compare '      super_obj{1,1}.label_signal{1,loop_over_assigned_signals} ' '    tmp_obj.signalf1{1,loop_over_peaks}])
                              %  disp(['compare '      num2str(size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},1)) ' '     num2str(size(tmp_obj.signalf1{1,loop_over_peaks},1))])
                              %  disp(['compare '       num2str(size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},2)) ' '     num2str(size(tmp_obj.signalf1{1,loop_over_peaks},2))])
                                if strcmp(super_obj{1,1}.label_signal{1,loop_over_assigned_signals} , tmp_obj.signalf1{1,loop_over_peaks})
                                    list_pos1=super_obj{1,1}.chemical_shift{1,loop_over_assigned_signals};
                                end
                            end
                            end
                        end
                         for loop_over_assigned_signals=1:size(super_obj{1,1}.label_signal,2)
                            if size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},2)>0
                            if size(tmp_obj.signalf2{1,loop_over_peaks},2)>0
                              %  disp(['compare '      super_obj{1,1}.label_signal{1,loop_over_assigned_signals} ' '    tmp_obj.signalf1{1,loop_over_peaks}])
                              %  disp(['compare '      num2str(size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},1)) ' '     num2str(size(tmp_obj.signalf1{1,loop_over_peaks},1))])
                              %  disp(['compare '       num2str(size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},2)) ' '     num2str(size(tmp_obj.signalf1{1,loop_over_peaks},2))])
                                if strcmp(super_obj{1,1}.label_signal{1,loop_over_assigned_signals} , tmp_obj.signalf2{1,loop_over_peaks})
                                    list_pos2=super_obj{1,1}.chemical_shift{1,loop_over_assigned_signals};
                                end
                            end
                            end
                         end
                         
                         if (list_pos2~=0) && (list_pos1~=0)
                             plot(list_pos2,list_pos1,[used_color '+'])
                             text(list_pos2,list_pos1,txt,'color',used_color)
                      min_chem_shift_spectrum1=min(min([min_chem_shift_spectrum1 list_pos1]));
                        max_chem_shift_spectrum1=max(max([max_chem_shift_spectrum1 list_pos1]));
        min_chem_shift_spectrum2=min(min([min_chem_shift_spectrum2 list_pos2]));
        max_chem_shift_spectrum2=max(max([max_chem_shift_spectrum2 list_pos2]));
                         end
    end
    margin=0.2;
    min_chem_shift_spectrum1=min_chem_shift_spectrum1 -margin;min_chem_shift_spectrum1=round(min_chem_shift_spectrum1*10)/10;
    max_chem_shift_spectrum1=max_chem_shift_spectrum1 +margin;max_chem_shift_spectrum1=round(max_chem_shift_spectrum1*10)/10;
    min_chem_shift_spectrum2=min_chem_shift_spectrum2 -margin;min_chem_shift_spectrum2=round(min_chem_shift_spectrum2*10)/10;
    max_chem_shift_spectrum2=max_chem_shift_spectrum2 +margin;max_chem_shift_spectrum2=round(max_chem_shift_spectrum2*10)/10;
    xlim([min_chem_shift_spectrum2 max_chem_shift_spectrum2])
    ylim([min_chem_shift_spectrum1 max_chem_shift_spectrum1])
    
    
    
    
    fig_num=list_2d(loop_over_spectra,1);
    drawnow
    if  exist('OCTAVE_VERSION', 'builtin') ~= 0
        print(['./O_spectrum_' dataset_name '_' tmp_obj.tag_name '.eps'],'-color');%for octave
    else
        print(['./M_spectrum_' dataset_name '_' tmp_obj.tag_name '.eps'],'-depsc');%for matlab
    end
end
end


