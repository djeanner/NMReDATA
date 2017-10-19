function [system,red_system]=get_from_dft_data(dft_folder,min_J,atoms_no_exchange,average_ch2,name)
average_ch3=1;%chemical shifts, copling, etc. for th three protons of CH3 will be averaged
if nargin<4
    average_ch2=0;
end
if nargin<3
    atoms_no_exchange=[6];
end
if nargin<2
    min_J=0.5;
end
if nargin<1
    path_many_hap='/Volumes/lacie_case/work_remove_duplicate_second/folder_1_234232341234434125467856/jeannerat_from_main_folder_oct_2015/Dropbox/HAP/gaussian_data/hap/';
    dft_folder=[path_many_hap '7-12dimethylbenzaanthracene.nmrJ.'];
end
%only for C-H (no OH !)
% table_1h_cs=zeros(0,0);
% table_coup_for_return=zeros(0,0);
disp_3d=1;
if disp_3d
    figure(888);clf;hold on;axis('equal');%pbaspect([1 1 1])
end
% dist_max_ch=1.6; %maximal distance of CH bond
% dist_max_cc=3; %maximal distance of CH bond

%% read dft files
system=read_dft_nmr_data_fn(dft_folder);
system.name=name;
%system has elements : atom_number, xyz, cs, J
%% create labels for atoms
for looop_of_x=1:size(system.atom_number,1)%for each C ...
    alias_table{looop_of_x}=num2str(looop_of_x);
    label_atoms{looop_of_x}=alias_table{looop_of_x};
end

system.bond_list=create_list_bonds(system.atom_number,system.xyz);
if disp_3d
    figure(889);clf;hold on;axis('equal');%pbaspect([1 1 1])
end
[cosy, hsqc, hmbc]=generate_correlations(system,atoms_no_exchange,disp_3d);
%% identify methyl
%% also make HSQC table and cosy 2J

%% add the number of the carbon to the label of the proton
number_of_h=zeros(1,size(system.atom_number,1));
for looop_of_x=1:size(system.atom_number,1)%for each C ...
    for heavy_atom=atoms_no_exchange
        if system.atom_number(looop_of_x,1)==heavy_atom%if carbon
            label_atoms{looop_of_x}=['(' alias_table{looop_of_x} ')'];
            atom_l='X';
            if heavy_atom==6 atom_l=['C'];  end
            if heavy_atom==7 atom_l=['N'];  end
            if heavy_atom==8 atom_l=['O'];  end
            for looop_of_h=1:size(system.atom_number,1)
                if system.atom_number(looop_of_h,1)==1% if H...
                    %                 tmp=xyz(looop_of_c,:)-xyz(looop_of_h,:);
                    %                 dist=sum(sum(tmp.*tmp));
                    %                 if dist<dist_max_ch
                    if system.bond_list(looop_of_h,looop_of_x)
                        number_of_h(  looop_of_x)=  number_of_h(  looop_of_x)+1;
                        label_atoms{looop_of_h}=['H' alias_table{looop_of_h} '(' atom_l  alias_table{looop_of_x} ')'];
                        
                    end
                end
            end
        end
    end
end


%% change lablels if file exists to give atom number labels (alias)
file=[ dft_folder 'log.alias'];
if exist(file,'file')
    fid=fopen(file);
    %    C = textscan(fid, '%s%s%f32%d8%u%f%f%s%f');
    C = textscan(fid, '%d%s');
    fclose(fid);
    numin=C{1};
    aliin=C{2};
    for loodpf=1:size(numin,1)
        tttpl=aliin(loodpf);
        label_atoms{numin(loodpf)}=tttpl{1};
    end
end

%% add labels to the 3D structure if was requested figure was initialized by generate_correlations
for looop_of_x=1:size(system.atom_number,1)%for each C ...
    if disp_3d
        text(system.xyz(looop_of_x,1) ,system.xyz(looop_of_x,2) ,system.xyz(looop_of_x,3) ,label_atoms(looop_of_x),'HorizontalAlignment','center','VerticalAlignment','middle')
    end
end
system.label_atoms=label_atoms;
clear label_atoms;

red_system=system;

%create lookup table to know in reduced system what atoms each signal
%correspodns to...
for loop_over_atoms=1:size(red_system.atom_number,1)
    red_system.lookup{loop_over_atoms}=[loop_over_atoms];
end
%% fuse protons of methyl and methylene if requested
OK=1;
while OK
    for looop_of_x=1:size(red_system.atom_number,1)
        if red_system.atom_number(looop_of_x,1)==6%if carbon
            count_h=0;list=[];
            for looop_of_h=1:size(red_system.atom_number,1)
                if red_system.atom_number(looop_of_h,1)==1% if H...
                    if red_system.bond_list(looop_of_h,looop_of_x)
                        count_h=count_h+1;
                        list=[list looop_of_h];
                    end
                end
            end
            if count_h==3
                if average_ch3
                    red_system=reduce_system(red_system,list);
                    red_system.label_atoms{list(1,1)}=['CH3' red_system.label_atoms{looop_of_x} ];
                    %                    disp(['removed protons for CH3:' num2str(looop_of_x)]);
                    break
                end
            end
            if count_h==2
                if average_ch2
                    red_system=reduce_system(red_system,list);
                    red_system.label_atoms{list(1,1)}=['CH2' red_system.label_atoms{looop_of_x} ];
                    %                    disp(['removed protons for CH2:' num2str(looop_of_x)]);
                    break;
                end
            end
        end
    end
    if looop_of_x>=size(red_system.atom_number,1) OK=0;end
end
if disp_3d
    figure(887);clf;hold on;axis('equal');%pbaspect([1 1 1])
end
[cosy, hsqc, hmbc]=generate_correlations(red_system,atoms_no_exchange,disp_3d);
%% add labels to the 3D structure if was requested figure was initialized by generate_correlations
for looop_of_x=1:size(red_system.atom_number,1)%for each C ...
    if disp_3d
        text(red_system.xyz(looop_of_x,1) ,red_system.xyz(looop_of_x,2) ,red_system.xyz(looop_of_x,3) ,red_system.label_atoms(looop_of_x),'HorizontalAlignment','center','VerticalAlignment','middle')
    end
end
number_of_different_h=sum(red_system.atom_number==1);




table_of_shift_ch3=zeros(1,size(red_system.atom_number,1));
table_1h_cs=zeros(1,number_of_different_h);
%table_1h_J=zeros(number_of_different_h,size(red_system.atom_number,1));
table_1h_J=zeros(number_of_different_h,number_of_different_h);
loop_over_1H=1;% in table of 1H only
%loop_over_all_1H=1;
text2='';
count_of_c=0;
table_ok=[];
for looop_of_x=1:size(red_system.atom_number,1)
    
    laelsh{looop_of_x}='';
end
for looop_of_x=1:size(red_system.atom_number,1)
    for heavy_atom=atoms_no_exchange
        if red_system.atom_number(looop_of_x,1)==heavy_atom
            count_table_of_shift_ch3=1;%table_of_shift_ch3
            count_of_c=count_of_c+1;
            textb=[ 'X'];
            if heavy_atom==6 textb=[ 'C'];  end
            if heavy_atom==7 textb=[ 'N'];  end
            if heavy_atom==8 textb=[ 'O'];  end
            here=number_of_h(looop_of_x);%1 for ch 2 for CH" 3 for CH3
            if here==0
                textb=[textb 'q '];
            end
            if here==1
                textb=[textb 'H '];
            end
            if here==2
                textb=[textb 'H2'];
            end
            if here==3
                textb=[textb 'H3'];
            end
            textb=[textb '(' num2str(count_of_c,'%3d') ')'];
            text2=[text2 textb ' '];
            
            for looop_of_h=1:size(red_system.atom_number,1)%for each C ...
                if red_system.atom_number(looop_of_h,1)==1
                    
                    %  tmp=xyz(looop_of_c,:)-xyz(looop_of_h,:);
                    %   dist=sum(sum(tmp.*tmp));
                    %   if dist<dist_max_ch,
                    if red_system.bond_list(looop_of_x,looop_of_h)
                        
                        
                        table_ok=[table_ok looop_of_h];
                        
                        table_1h_cs(1,loop_over_1H)=red_system.cs(looop_of_h);
                        text2=[text2 num2str(table_1h_cs(1,loop_over_1H),': %.2f') ','];
                        laelsh{loop_over_1H}=[textb num2str(red_system.cs(looop_of_x,1),'%.3f') ' '];
                        
                        num_h=1;
                        for looop_of_h2=1:size(red_system.atom_number,1)%for each C ...
                            % if atom_number(looop_of_h2,1)==1,
                            table_1h_J(loop_over_1H , num_h)=red_system.J(looop_of_h,looop_of_h2);
                            num_h=num_h+1;
                            % end
                        end
                        ref_to_atom(loop_over_1H,1)=looop_of_h;
                        
                        loop_over_1H=loop_over_1H+1;
                        %  loop_over_all_1H=loop_over_all_1H+1;
                    end
                    
                end
            end
            
        end
        text2=[text2 '\n'];
    end
end
%sprintf(text2)

table_ok=1:number_of_different_h;


%table_of_shift_ch3;
table_ok2=[];
for loo1=1:size(red_system.atom_number,1)%  number_of_different_h
    if table_of_shift_ch3(loo1)==0
        if red_system.atom_number(loo1,1)==1
            table_ok2=[table_ok2 loo1];
        end
    end
end
for loo1=1:number_of_different_h%  number_of_different_h
    
    %  for loo2=1:size(atom_number,1),%  number_of_different_h
    %  if   abs(table_1h_J(loo1,loo2))>min_J,
    % laelsh{loo1}= [ laelsh{loo1} ' ' num2str(table_1h_J(loo1,loo2),'%.1f') ' '];
    %    if atom_number(loo2,1)==1,
    %   laelsh{loo1}= [ laelsh{loo1}  '(' num2str(loo1) ',' num2str(loo2) ')' num2str(table_1h_J(loo1,loo2),'%.1f') ' '];
    
    %     text2=[text2 '(' num2str(loo1) ',' num2str(loo2) ')' num2str(table_1h_J(loo1,loo2),'%.1f') ' '];
    % end
    % end
    % end
    %for loo1=1:size(atom_number,1),%
    if table_of_shift_ch3(loo1)~=0,
        work_on=table_of_shift_ch3(loo1);
        for loo2=1:size(red_system.atom_number,1),%
            if work_on==table_of_shift_ch3(loo2),
                for loo3=1:number_of_different_h,
                    table_1h_J(loo3,work_on)=table_1h_J(loo3,work_on)+table_1h_J(loo3,loo2);
                end
                table_of_shift_ch3(loo2)=0;
            end
            
        end
        for loo3=1:number_of_different_h,
            table_1h_J(loo3,work_on)=table_1h_J(loo3,work_on)/3;
            %  disppp=table_1h_J(loo3,work_on)
        end
    end
end

table_1h_J=table_1h_J(:,table_ok);




%text
text2='';
table_coup_for_return=zeros(number_of_different_h,number_of_different_h);
for loo1=1:size(table_1h_J,1)
    %for loo1=1:size(laelsh)
    inc2=1;
    text2=[text2 laelsh{loo1}] ;
    text2=[text2 num2str(table_1h_cs(1,loo1),'%.3f') ' J= '];
    
    for loo2=1:size(table_1h_J,2)
        if abs(table_1h_J(loo1,loo2))>min_J
            % text2=[text2 num2str(table_1h_J(loo1,loo2),'%.1f') ' '];
            text2=[text2 '(' num2str(loo1) ',' num2str(loo2) ') ' num2str(table_1h_J(loo1,loo2),'%.1f') ' '];
            %    text2=[text2 '(' num2str(loo2) ',' num2str(loo1) ')' num2str(table_1h_J(loo2,loo1),'%.1f') ' '];
            table_coup_for_return(loo1,inc2)=table_1h_J(loo1,loo2);
            inc2=inc2+1;
        end
    end
    
    text2=[text2 '\n'];
    
end
%
% end
read_dft_nmr_data_fn_disp= sprintf(text2);
disp(read_dft_nmr_data_fn_disp);


table_coup_for_return=abs(table_coup_for_return);

if disp_3d
    if exist('savefig','builtin')
        savefig([ 'Structure_with_COSY_HSQC_HMBC_correlations_' name '.fig'])
    end
end
label_atoms=red_system.label_atoms;
nmrJ_log=red_system.J;
red_system.nb_bond_between=create_nb_bonds_between_atoms(red_system.bond_list);
save_sdf(system,red_system,[system.name '.sdf']);
save_sdf(system,red_system,[system.name '_implicit.sdf'],'implicit');


end
