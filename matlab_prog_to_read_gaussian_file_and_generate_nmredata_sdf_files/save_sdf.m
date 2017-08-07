function save_sdf(system,red_system,file_name,varargin)
field_delimiter=', ';%', ' or only ','
correl_delimiter='/';%', ' or only ','
factor_change_unit_for_coordinates=1;
dim_order=[1:3];%for xyz corrdinates used to project differently
pretag_label='NMREDATA_';
larmorh=400;%default value
larmorc=larmorh*67.262/267.513;%default value
larmorn=larmorh*19.331/267.513;%default value
list_hh=1;%to list J(H,H) in tag J
list_ch=1;%to list J(C,H) in tag J
list_cc=0;%to list J(C,C) in tag J
list_all=0;%to list all J(C,C) in tag J

assign_jhh=1;%to list partnerrs of couling in 1D 1H


comment=1;% this is to dispay or not comments
minj=1;% defaults minimal coupling considered value may be changed (see below)
start_of_comment_char=';';% this is the the caracte starting comment
chem_el='H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHzTlPbBiPoAtRnFrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDsRg';
explit_h=1;
for i = 1:length(varargin)
    cur_arg = varargin{i};
    if strcmp(cur_arg,'list_CC') list_cc=1;
    elseif strcmp(cur_arg,'implicit') explit_h=0;
    elseif strcmp(cur_arg,'list_CH') list_ch=1;
    elseif strcmp(cur_arg,'list_HH') list_hh=1;
    elseif strcmp(cur_arg,'nocomment') comment=0;
    elseif strcmp(cur_arg,'minJ') minj=varargin{i+1};
    elseif strcmp(cur_arg,'larmorh') larmorh=varargin{i+1};larmorc=larmorh*67.262/267.513; larmorn=larmorh*19.331/267.513;%default value
    elseif strcmp(cur_arg,'larmorc') larmorch=varargin{i+1};
    elseif strcmp(cur_arg,'larmorn') larmornh=varargin{i+1};
    end
end

%% save sdf file version 0.93
if nargin==0
    tmp='./for_test_system.mat';
    load(tmp)
    file_name='./del.sdf';
else
    tmp='./for_test_system.mat';
    %save(tmp,'system')
    save(tmp)
end

% prepare molecule part of .sdf file
line2='Demo .sdf files generated from gaussian DFT/GIAO by matlab functions written by Damien Jeannerat in April 2017';%origin of the .sdf file
line3='dj_test_generationV_0.93 1';
nb_labels=size(red_system.atom_number,1);

nb_atoms=size(system.atom_number,1);
if explit_h
    list_of_atoms_in_mol_file=1:nb_atoms;% all the atoms
else
    nb_atoms=sum(system.atom_number~=1);%number of non-H atoms
    list_of_atoms_in_mol_file=find(system.atom_number~=1)';% non-H atoms
end
inc=1;
pointer=zeros(1,max(max(list_of_atoms_in_mol_file)));
for loop_over_atoms=list_of_atoms_in_mol_file
    pointer(1,loop_over_atoms)=inc;inc=inc+1;
end
nb_bond=sum(sum(system.bond_list(list_of_atoms_in_mol_file,list_of_atoms_in_mol_file)>0))/2;
%to be verfied...
%el=1;disp(['Element number ' num2str(el) ':' chem_el(el*2-1)  chem_el(el*2)]);% demo use chem_el
% el=1;txt1=chem_el(el*2-1);  if chem_el(el*2)~=32 txt1=([txt1 chem_el(el*2)]);end%

fid=fopen(file_name,'w');
if fid==-1%OK
    disp(['Could not write file: ' file_name] );
else
    
    write_mole_part_of_file;% not working well!!
    write_tag_first_few_tags;%this function is embedded to this code (see below)
    write_tag_nmr_signals;%this function is embedded to this code (see below)
    write_tag_j;%this function is embedded to this code (see below)
    write_nmr_data_1D(1);%this function is embedded to this code (see below)
    write_nmr_data_1D(6);%this function is embedded to this code (see below)
    %these functions are embedded to this code (see below)
    write_nmr_data_1D(6,'dept135');
    write_nmr_data_2D(1,1,'N','J',1);%cosy if add a coupling, will take all coupling larger than the value otherwise only 2J adn 3J
    write_nmr_data_2D(6,1,'1','J');%hsqc
    write_nmr_data_2D(6,1,'N','J',1);%hmbc if add a coupling, will take all coupling larger than the value otherwise only 2J adn 3J
    write_origin;
    four_dol='$$$$';
    fprintf(fid,'%s\n',four_dol);%closing record (there may be more than on structure in .sdf files)
    fclose(fid);
end
    function write_nmr_data_2D(element1,element2,range_mixing,mixing,min_val)
        if nargin<5
            min_val=0;
        end
        
        iso1_txt='';
        iso2_txt='';
        if element1==1,iso1_txt='1H';larmor=larmorh;end
        if element1==6,iso1_txt='13C';larmor=larmorc;end
        if element1==7,iso1_txt='15N';larmor=larmorn;end
        if element2==1,iso2_txt='1H';larmor=larmorh;end
        if element2==6,iso2_txt='13C';larmor=larmorc;end
        if element2==7,iso2_txt='15N';larmor=larmorn;end
        
        fprintf(fid,['>  <' pretag_label '2D_' iso1_txt '_' range_mixing  mixing '_' iso2_txt '>\n']);
        if min_val==0
            fprintf(fid,';list correlations when the number of bonds is 1 for HSQC, 2-3 for COSY and HMBC\n');
        else
            fprintf(fid,';list correlations when abs(J) > %.2f \n',min_val);
        end
        fprintf(fid,'Larmor=%f\n',larmor);
        cortype='';
        if (strcmp(range_mixing,'1') && strcmp(mixing,'J') && (((element1==6) && (element2==1))||((element1==1) && (element2==6)))) cortype='HSQC';end
        if (strcmp(range_mixing,'N') && strcmp(mixing,'J') && (((element1==6) && (element2==1))||((element1==1) && (element2==6)))) cortype='HMBC';end
        if (strcmp(range_mixing,'N') && strcmp(mixing,'J') && (element1==1) && (element2==1)) cortype='COSY';end
        fprintf(fid,['CorType=' cortype '\n']);
        range_values=[];
        if strcmp(range_mixing,'1');range_values=[1];end
        if strcmp(range_mixing,'N');range_values=[2 3];end
        for loop_over_atoms1=1:nb_labels
            el1=red_system.atom_number(loop_over_atoms1,1);
            if el1==element1
                for loop_over_atoms2=(element1==element2)*loop_over_atoms1+1:nb_labels
                    el2=red_system.atom_number(loop_over_atoms2,1);
                    if el2==element2
                        if loop_over_atoms1~=loop_over_atoms2
                            if min_val==0
                                crit=(size(find(range_values==red_system.nb_bond_between(loop_over_atoms1,loop_over_atoms2)),2)>0);
                            else
                                crit=abs(red_system.J(loop_over_atoms1,loop_over_atoms2))>min_val;
                            end
                            if crit
                                % list the correlation (mandatory)
                                fprintf(fid,'%s%s%s',red_system.label_atoms{loop_over_atoms1},correl_delimiter,red_system.label_atoms{loop_over_atoms2});
                                % active coupling (coupling betweent
                                % the correlated spins)
                                fprintf(fid,'%sJa=%.2f',field_delimiter,red_system.J(loop_over_atoms1,loop_over_atoms2));
                                % F1 passive coupling (coupling the spin in F1 and its partners)
                                first=1;
                                for other=1:nb_labels
                                    if other~=loop_over_atoms1
                                        if other~=loop_over_atoms2
                                            if red_system.atom_number(loop_over_atoms1,1)==1
                                                if red_system.atom_number(other,1)==1
                                                    crit=abs(red_system.J(loop_over_atoms1,other))>minj;
                                                    if crit
                                                        if first
                                                        fprintf(fid,'%sJ1=%.2f(%s)',field_delimiter,red_system.J(loop_over_atoms1,other),red_system.label_atoms{other});   first=0;

                                                        else
                                                            
                                                        fprintf(fid,'%s%.2f(%s)',field_delimiter,red_system.J(loop_over_atoms1,other),red_system.label_atoms{other});
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                                first=1;
                                
                                % F2 passive coupling (coupling the spin in F1 and its partners)
                                for other=1:nb_labels
                                    if other~=loop_over_atoms1
                                        if other~=loop_over_atoms2
                                            if red_system.atom_number(loop_over_atoms2,1)==1
                                                if red_system.atom_number(other,1)==1
                                                    crit=abs(red_system.J(loop_over_atoms2,other))>minj;
                                                    if crit
                                                        %  fprintf(fid,'%sJ2(%s)=%.2f',field_delimiter,red_system.label_atoms{other},red_system.J(loop_over_atoms2,other));
                                                        
                                                        if first
                                                            fprintf(fid,'%sJ2=%.2f(%s)',field_delimiter,red_system.J(loop_over_atoms2,other),red_system.label_atoms{other});first=0;
                                                        else
                                                            
                                                            fprintf(fid,'%s%.2f(%s)',field_delimiter,red_system.J(loop_over_atoms2,other),red_system.label_atoms{other});
                                                            
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                                if comment
                                    
                                    fprintf(fid,start_of_comment_char);% comment starts here
                                    if isfield(red_system,'nb_bond_between')% this will all allow to specify 1J, 2J, etc.
                                        fprintf(fid,'%d',red_system.nb_bond_between(loop_over_atoms1,loop_over_atoms2));
                                    end
                                    %labeling the element
                                    txt1=chem_el(el1*2-1);  if chem_el(el1*2)~=32 txt1=([txt1 chem_el(el1*2)]);end%
                                    txt2=chem_el(el2*2-1);  if chem_el(el2*2)~=32 txt2=([txt2 chem_el(el2*2)]);end%
                                    fprintf(fid,'J(%s,%s)',txt1,txt2);
                                    if isfield(red_system,'nb_bond_between')% this will all allow to specify 1J, 2J, etc.
                                        
                                        if red_system.nb_bond_between(loop_over_atoms1,loop_over_atoms2)>3
                                            fprintf(fid,' Note this is a >3J ');
                                        end
                                    end
                                    
                                end
                                fprintf(fid,'\n');
                                
                            end
                        end
                    end
                end
                
                
                %                 if length(mult_struc_text)==0 mult_struc_text='s';end
                %                 fprintf(fid,'%sS=%s',field_delimiter,mult_struc_text);% chemical shift and label
                %                 if length(list_coupling_text)~=0
                %                     fprintf(fid,'%sJ=%s',field_delimiter,list_coupling_text);% chemical shift and label
                %                 end
                %                 fprintf(fid,'%sE=%d',field_delimiter,size(red_system.lookup{loop_over_atoms1},2));% chemical shift and label
                %                 fprintf(fid,'\n');
            end
        end
        fprintf(fid,'\n');%end of tag
    end
    function write_nmr_data_1D(element,type)
        if nargin<2
            type='';
        end
        
        iso_txt='';
        dec_1h=0;
        if element==1,iso_txt='1H';larmor=larmorh;end
        if element==6,iso_txt='13C';larmor=larmorc;dec_1h=1;end
        if element==7,iso_txt='15N';larmor=larmorn;dec_1h=1;end
        fprintf(fid,['>  <' pretag_label '1D_' iso_txt '>\n']);
        if element==1
            fprintf(fid,[';J listed when the H-H couplings >' num2str(minj) 'Hz\n']);
        end
        fprintf(fid,'Larmor=%f\n',larmor);
        if strcmp(type,'dept135') dept135=1;else dept135=0;end
        if dec_1h
            fprintf(fid,'Decoupled=1H\n');
        end
        if dept135
            fprintf(fid,'Sequence=DEPT135\n');
        end
        inc=1;table_to_sort_list=zeros(1,nb_labels);
        for loop_over_atoms1=1:nb_labels
            el1=red_system.atom_number(loop_over_atoms1,1);
            if el1==element
                if el1==1%for 1H spectra list JHH and give multiplicity
                    % fprintf(fid,'%.2f%sL=%s',red_system.cs(loop_over_atoms1,1),field_delimiter,red_system.label_atoms{loop_over_atoms1});% chemical shift and label
                    tmp1=sprintf('%.4f%sL=%s',red_system.cs(loop_over_atoms1,1),field_delimiter,red_system.label_atoms{loop_over_atoms1});% chemical shift and label
                    table_to_sort_list(1,inc)=red_system.cs(loop_over_atoms1,1);
                    %study partners
                    mult_struc_text='';
                    list_coupling_text='';                        count=0;
                    
                    for loop_over_atoms2=1:nb_labels
                        if loop_over_atoms1~=loop_over_atoms2
                            if red_system.atom_number(loop_over_atoms2,1)==1
                                if abs(red_system.J(loop_over_atoms1,loop_over_atoms2))>minj% is coupling large enough to be considered
                                    
                                    toto=red_system.lookup{loop_over_atoms2};
                                    if size(toto,2)==1 lab='d';end
                                    if size(toto,2)==2 lab='t';end
                                    if size(toto,2)==3 lab='q';end
                                    if size(toto,2)>3 lab='DEFINE_ME_HERE';end
                                    mult_struc_text=[mult_struc_text lab];
                                    if count
                                        list_coupling_text=[list_coupling_text field_delimiter ];
                                    end
                                    list_coupling_text=[list_coupling_text  num2str(red_system.J(loop_over_atoms1,loop_over_atoms2),'%.2f')];
                                    if assign_jhh
                                        list_coupling_text=[list_coupling_text '(' red_system.label_atoms{loop_over_atoms2} ')'];
                                    end
                                    count=count+1;
                                end
                            end
                        end
                    end
                    if length(mult_struc_text)==0 mult_struc_text='s';end
                    tmp2=sprintf('%sS=%s',field_delimiter,mult_struc_text);% chemical shift and label
                    tmp3='';
                    if length(list_coupling_text)~=0
                        tmp3=sprintf('%sJ=%s',field_delimiter,list_coupling_text);% chemical shift and label
                    end
                    tmp4=sprintf('%sE=%d',field_delimiter,size(red_system.lookup{loop_over_atoms1},2));% chemical shift and label
                    text_for_line{inc}=[tmp1 tmp2 tmp3 tmp4];
                    inc=inc+1;
                    %fprintf(fid,'%s\n',text_for_line{inc});inc=inc+1;
                    
                end
                if el1==6%for 13H
                    %determine CH CH2 CH3 for intensity
                    nb_h=0;
                    for loop_over_atoms2=1:nb_labels
                        if loop_over_atoms1~=loop_over_atoms2
                            if red_system.bond_list(loop_over_atoms1,loop_over_atoms2)
                                if red_system.atom_number(loop_over_atoms2,1)==1
                                    nb_h=nb_h+size(red_system.lookup{loop_over_atoms2},2);
                                end
                            end
                        end
                    end
                    if (dept135 && (nb_h==2)) sign_of_signal=-1;else sign_of_signal=1;end
                    if (dept135 && (nb_h==0)) sign_of_signal=0;end
                    if sign_of_signal~=0,
                        tmp1=sprintf('%.4f%sL=%s',red_system.cs(loop_over_atoms1,1),field_delimiter,red_system.label_atoms{loop_over_atoms1});% chemical shift and label
                        table_to_sort_list(1,inc)=red_system.cs(loop_over_atoms1,1);
                        
                        tmp2=sprintf('%sI=%.2f',field_delimiter,sign_of_signal*(100+10*nb_h));% chemical shift and label
                        text_for_line{inc}=[tmp1 tmp2 ];
                        inc=inc+1;
                        
                    end
                end
                
                
            end
        end
        [values out_list]=sort(table_to_sort_list(1,1:inc-1),'descend');
        
        for list_order=1:inc-1
            fprintf(fid,'%s\n',text_for_line{out_list(list_order)});
        end
        fprintf(fid,'\n');%end of tag
    end


    function write_tag_nmr_signals
        
        fprintf(fid,['>  <' pretag_label 'SIGNALS>\n']);
        for loop_over_atoms=1:nb_labels
            fprintf(fid,'%s%s%.4f',red_system.label_atoms{loop_over_atoms},field_delimiter,red_system.cs(loop_over_atoms,1));
            for loop_over_list_of_atoms=red_system.lookup{loop_over_atoms}
                if explit_h
                    fprintf(fid,'%s%d',field_delimiter,loop_over_list_of_atoms);
                else
                    if red_system.atom_number(loop_over_atoms)~=1
                        fprintf(fid,'%s%d',field_delimiter,pointer(1,loop_over_list_of_atoms));
                    else
                        refc=1;
                        for loop_over_atoms2=1:nb_labels
                            if red_system.bond_list(loop_over_atoms,loop_over_atoms2)
                                refc=loop_over_atoms2;
                            end
                        end
                        refc=red_system.lookup{refc};
                        refc=refc(1,1);
                        fprintf(fid,'%sH%d',field_delimiter,pointer(1,refc));
                        break;
                        
                    end
                end
            end
            fprintf(fid,'\n');
            
        end
        fprintf(fid,'\n');%end of tag
        
    end
    function write_tag_j
        
        %% write table of coupling constants (full coupling network
        fprintf(fid,['>  <' pretag_label 'J>\n']);
        %comment about what is listed
        if list_all
            fprintf(fid,[';list all couplings >' num2str(minj) 'Hz\n']);
        else
            if list_hh fprintf(fid,[';list abs(J(H,H))>' num2str(minj) 'Hz\n']);end% comment
            if list_ch fprintf(fid,[';list abs(J(C,H))>' num2str(minj) 'Hz\n']);end% comment
            if list_cc fprintf(fid,[';list abs(J(C,C))>' num2str(minj) 'Hz\n']);end% comment
        end
        
        for loop_over_atoms1=1:nb_labels
            for loop_over_atoms2=loop_over_atoms1+1:nb_labels
                if abs(red_system.J(loop_over_atoms1,loop_over_atoms2))>minj% is coupling large enough to be considered
                    should_list_this_coupling=list_all;% to only list some types of couplings
                    el1=red_system.atom_number(loop_over_atoms1,1);
                    el2=red_system.atom_number(loop_over_atoms2,1);
                    if el2>el1
                        tmp_val=el1;el1=el2;el2=tmp_val;
                        s2=loop_over_atoms1;
                        s1=loop_over_atoms2;
                    else
                        s1=loop_over_atoms1;
                        s2=loop_over_atoms2;
                    end
                    if (list_hh && (el1==1) && (el2==1)) should_list_this_coupling=1;end
                    if (list_cc && (el1==6) && (el2==6)) should_list_this_coupling=1;end
                    if (list_ch && (el1==6) && (el2==1)) should_list_this_coupling=1;end
                    if should_list_this_coupling
                        fprintf(fid,'%s%s%s%s%.2f',red_system.label_atoms{s1},field_delimiter,red_system.label_atoms{s2},field_delimiter,red_system.J(loop_over_atoms1,loop_over_atoms2));
                        % comment specify the type of coupling
                        % nJ(isotope1,isotope2)
                        if comment % add the type of coupling 2JCH, ...
                            fprintf(fid,start_of_comment_char);% comment starts here
                            if isfield(red_system,'nb_bond_between')% this will all allow to specify 1J, 2J, etc.
                                fprintf(fid,'%d',red_system.nb_bond_between(loop_over_atoms1,loop_over_atoms2));
                            end
                            %labeling the element
                            txt1=chem_el(el1*2-1);  if chem_el(el1*2)~=32 txt1=([txt1 chem_el(el1*2)]);end%
                            txt2=chem_el(el2*2-1);  if chem_el(el2*2)~=32 txt2=([txt2 chem_el(el2*2)]);end%
                            fprintf(fid,'J(%s,%s)',txt1,txt2);
                            if isfield(red_system,'nb_bond_between')% this will all allow to specify 1J, 2J, etc.
                                
                                if red_system.nb_bond_between(loop_over_atoms1,loop_over_atoms2)>3
                                    fprintf(fid,' Note this is a >3J ');
                                end
                            end
                        end
                        fprintf(fid,'\n');
                    end
                end
            end
        end
        fprintf(fid,'\n');%end of tag
        
    end
    function write_tag_first_few_tags
        fprintf(fid,['>  <' pretag_label 'VERSION>\n']);
        fprintf(fid,'0.93\n');
        fprintf(fid,'\n');%end of tag
        
        fprintf(fid,['>  <' pretag_label 'ID>\n']);
        fprintf(fid,'no ide yet\n');
        fprintf(fid,'\n');%end of tag
        
        fprintf(fid,['>  <' pretag_label 'SOLVENT>\n']);
        fprintf(fid,'cdcl3\n');
        fprintf(fid,'\n');%end of tag
        
        fprintf(fid,['>  <' pretag_label 'TEMPERATURE>\n']);
        fprintf(fid,'273.15\n');
        fprintf(fid,'\n');%end of tag
        
    end
    function write_mole_part_of_file
        fprintf(fid,'%s\n',system.name);
        fprintf(fid,'%s\n',line2);
        fprintf(fid,'%s\n',[line3 ' ' date ]);
        fprintf(fid,'%2d %2d  0  0  0  0  0  0  0  0999 V2000\n',nb_atoms,nb_bond);
        %% problems : the mol part is not OK. No double bond (important for chem draw not to transform C=O in to C-O-H...
        %% problems : the mol part is not OK. Problem with scales
        inc=1;pointer=zeros(1,max(max(list_of_atoms_in_mol_file)));
        for loop_over_atoms=list_of_atoms_in_mol_file
            el=system.atom_number(loop_over_atoms,1);
            fprintf(fid,' %9.4f %9.4f %9.4f %s%s  0  0  0  0  0  0  0  0  0  0  0  0\n',factor_change_unit_for_coordinates*system.xyz(loop_over_atoms,dim_order(1,1)),factor_change_unit_for_coordinates*system.xyz(loop_over_atoms,dim_order(1,2)),factor_change_unit_for_coordinates*system.xyz(loop_over_atoms,dim_order(1,3)),chem_el(el*2-1),chem_el(el*2));
            pointer(loop_over_atoms)=inc;inc=inc+1;
        end
        for loop_over_atoms=list_of_atoms_in_mol_file
            for loop_over_atoms2=list_of_atoms_in_mol_file
                if loop_over_atoms<loop_over_atoms2
                    nb=system.bond_list(loop_over_atoms,loop_over_atoms2);
                    if nb
                        fprintf(fid,'%2d %2d %2d  0  0  0\n',pointer(loop_over_atoms),pointer(loop_over_atoms2),nb);
                    end
                end
            end
        end
        fprintf(fid,'M  END\n');%end of molecule part 1/2
        fprintf(fid,'\n');%end of molecule part 2/2
    end
    function write_origin
        if isfield(red_system,'origin')
            fprintf(fid,['>  <' pretag_label 'ORIGIN>\n']);
            if contains(red_system.origin,'dft_simulation')
                fprintf(fid,['Source=Calculation\n']);
                fprintf(fid,['Geometry=%s\n'],red_system.dft_param_geometry);
                fprintf(fid,['Shielding=%s\n'],red_system.dft_param_nmrcs);
                fprintf(fid,['Coupling=%s\n'],red_system.dft_param_nmrj');
                fprintf(fid,['Software=\n'],red_system.dft_param_software);
                fprintf(fid,['Version=%s\n'],red_system.dft_param_version);
            end
            fprintf(fid,'\n');%end of tag
        end
    end
% 2,2'-[(2,6-Di-tert-butyl)anthracen-9,10-diyl]dithiophen
%   CDK     0227171006
% nmrshiftdb2 30001243
%  32 36  0  0  0  0  0  0  0  0999 V2000
%    -2.8438    0.9188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0

%...
%   1  2  1  0  0  0  0
%...
%  32 30  1  0  0  0  0
% M  END
%
% >  <VERSION>
% 1.0
%
% >  <ID>
% nmrshiftdb2:30001243
%
%
%
% >  <ORIGIN>
% Source=Literature
% Reference=no doi exists
% CompoundNumber=41


end