function [spinach_system, optim_param, simple_data, toto_check]=generate_spinach_spin_system(list_of_nmredata, item_num, iso, list_methyl, list_methylene, min_j, keep_label, keep_label_core )
% should read data from spectra, not from assignemnt and J tags... (for the
% case when they are not assigned...
if nargin<6
    min_j=1;
end

verbose=3;
list_diff_J1=[];
list_diff_J2=[];
simple_data.main_cs_table=[];
simple_data.main_cs_table_lin=[];
simple_data.main_j_table2d=[];
bas=struct;
counter_number_sym_group_magn_equ_spin=1;
% bas.sym_group={'S3','S3','S3'};
% bas.sym_spins={[14 15 16],[17 18 19],[20 21 22]};

if nargin<3
    iso='1H';% element for which atoms spin system is generated
end
if nargin<4
    list_methyl=[];
end
if nargin<5
    list_methylene=[];
end
if nargin<7
    keep_label=ones(1,10000);
    %     else
    %         keep_label=list_keep_label{looop_over_compounds};
end

if nargin<8
    keep_label_core=ones(1,10000);
    %     else
    %         keep_label_core=list_keep_label_core{looop_over_compounds};
end

%convert iso into element
for i=1:size(iso,2)
    if iso(1,i)
        tmp=strfind('0123456789',iso(1,i));
        islet(1,i)=1;
        
        if size(tmp,2)>0
            if tmp(1,1)>0
                islet(1,i)=0;
                
            end
        end
    end
end
element=iso(find(islet==1));

inter.zeeman.scalar={};
inter.coupling.scalar={};
%  inter.coupling=struct();
%inter.coupling=struct();

lookup_cs_compound=[];
lookup_cs_atom=[];
lookup_j_compound=[];
lookup_j_atom=[];
list_diff_chem_shift_returned =[];
filled_by_other_compounds_chemical_shifts=0;
filled_by_other_compounds_j=0;
filled_by_other_compounds_cs=0;
filled_by_other_compounds_cs2=0;
filled_large=0;
nb_molecules=size(item_num,2);
start_pos_compound_cs=zeros(1,nb_molecules);
for looop_over_compounds=1:nb_molecules
    if verbose==3
        disp(['Compound : ' num2str(looop_over_compounds) ' out of ' num2str(nb_molecules) '========================================' ])
    end
    if item_num(1,looop_over_compounds)>0
        nmredata=list_of_nmredata{looop_over_compounds};
        tmi=item_num(1,looop_over_compounds);
        current_obj=nmredata{tmi};
        larmor=current_obj.larmor;%should be the same... but never know could mix spectrometers
        counter_param_coup=1;%couter couplings.
        counter_param_spin_returned=1;%couter chem shifts
        
        counter_spins=1;%% MM
        
        
        list_diff_chem_shift=[];counter_param_spin=1;
        
        counter_all_lables_including_skiped_ones=1;
        
        %tot.dfs=struct()
        %  inter.coupling.scalar{ll1,ll2}=0;
        
        
        % % convert element into isotope
        % iso='';
        % if strcmp('H',element)%test if is hydrogen
        %     iso='1H';
        % end
        % if strcmp('C',element)%test if is carbon
        %     iso='13C';
        %     %  iso='1H';
        % end
        % if strcmp('P',element)%test if is P
        %     iso='31P';
        % end
        % if strcmp('F',element)%test if is F
        %     iso='19F';
        % end
        % if size(iso,2)==0
        %     error('isotope not identified')
        % end
        
        %% get assignment tag
        [assignment_obj,j_obj]=get_assignement_objects(nmredata);
        
        
        if verbose
            disp('Only H with labels and chemical shift in the NMREDATA_ASSIGNMENT are considered')
        end
        
        %% start generation of spinach data
        
        %determine larmor frequency
        % sys.magnet=2*pi*500.101412e6/spin('1H');% "spi"n a spinach function in Tesla
        sys.magnet=2*pi*larmor*1e6/spin(iso);% "spi"n a spinach function in Tesla
        if verbose
            disp(['Magnet : ' num2str(sys.magnet) ' Tesla (used by spinach) Calculation based on Larmor=' num2str(larmor) ' MHz for isotope : ' iso  ])
        end
        for lo=1:size(assignment_obj.atom_num,1)% loop over atoms
            cur_atom_num=assignment_obj.atom_num{lo,1};
            this_is_correct_element=0;
            how_many_atoms=1;%may be increased if implicit H
            if cur_atom_num(1,1)<0% explic H to this atom....
                if strcmp('H',element)%test if is hydrogen
                    this_is_correct_element=1;
                    if verbose>2
                        disp('found implicit H');
                        %   how_many_hydrogen_atom_are_linked_to_this_element=determine_how_many_hydrogen_atom_are_linked_to_this_element(assignment_obj.structure,cur_atom_num(1,1));
                    end
                    warning(['how_many_hydrogen_atom_are_linked_to_this_element is it CH3, CH2, CH ' num2str(cur_atom_num) ])
                    how_many_atoms=1;
                    
                    %list methyl
                    if size(find(list_methyl==-cur_atom_num),2)
                        how_many_atoms=3;
                    end
                    
                    %list methylene
                    if size(find(list_methylene==-cur_atom_num),2)
                        how_many_atoms=2;
                    end
                    
                    disp(['set to ' num2str(how_many_atoms)])
                end
            else
                atom_cell=assignment_obj.structure.atom.n(1,cur_atom_num(1,1));
                at_str=atom_cell{1,1};
                
                if strcmp(at_str,element)
                    this_is_correct_element=1;
                    if verbose>2
                        disp('found atom');
                    end
                end
            end
            
            if keep_label(1,counter_all_lables_including_skiped_ones)%only if in list of kept labels
                if this_is_correct_element
                    
                    loc_list_atom_number=[];
                    list_label_number(counter_param_spin)=    lo;
                    labcur=assignment_obj.label_signal(1,lo);
                    for loop=1:how_many_atoms% usually 1, but 2 for CH2 and 3 for CH3
                        sys.isotopes{1,counter_spins+filled_by_other_compounds_cs}=iso;%HERE may not be iso... but partner of J...
                        if size(item_num,2)>1
                            label=[num2str(looop_over_compounds) ':' labcur{1,1}];
                        else
                            label=[ labcur{1,1}];
                        end
                        sys.labels{1,counter_spins+filled_by_other_compounds_cs}=label;
                        inter.zeeman.scalar{1,counter_spins+filled_by_other_compounds_cs}=assignment_obj.chemical_shift{1,lo};
                        %    no=get_list_of_spin_number(assignment_obj,cur_atom_lab);
                        list_diff_chem_shift(counter_param_spin,loop)=counter_spins;
                        
                        if keep_label_core(1,counter_all_lables_including_skiped_ones)
                            
                            list_diff_chem_shift_returned(filled_by_other_compounds_chemical_shifts+counter_param_spin_returned,loop)=counter_spins;
                            lookup_cs_atom(filled_by_other_compounds_chemical_shifts+counter_param_spin_returned)=counter_param_spin_returned;
                            lookup_cs_compound(filled_by_other_compounds_chemical_shifts+counter_param_spin_returned)=looop_over_compounds;
                            
                            
                        end
                        loc_list_atom_number=[loc_list_atom_number counter_spins];
                       %%  simple_data.main_cs_table(1,filled_large+lo)=assignment_obj.chemical_shift{1,lo};%
            simple_data.main_cs_table(  looop_over_compounds,lo)=assignment_obj.chemical_shift{1,lo};%&&&&&&&&&&&&&&&&&&
                        simple_data.main_cs_table_lin(  1,filled_by_other_compounds_cs+counter_spins)=assignment_obj.chemical_shift{1,lo};%&&&&&&&&&&&&&&&&&&
            simple_data.main_cs_table2(  looop_over_compounds,counter_spins)=assignment_obj.chemical_shift{1,lo};%
                        simple_data.lookup_cs_table(looop_over_compounds,0*filled_by_other_compounds_cs+counter_spins)=lo;
                        simple_data.lookup_cs_table_lin                 (1,filled_by_other_compounds_cs+counter_spins)=lo;
                         simple_data.lookup_cs_table2(looop_over_compounds,0*filled_by_other_compounds_cs+counter_spins)=counter_spins;
                         simple_data.lookup_cs_table_lin2                 (1,filled_by_other_compounds_cs+counter_spins)=counter_spins;
                        simple_data.lookup_compound_table               (1,filled_by_other_compounds_cs+counter_spins)=looop_over_compounds;
                        simple_data.label{                                 filled_by_other_compounds_cs+counter_spins}=label;
                  
% simple_data.main_cs_table2(  looop_over_compounds,lo)=assignment_obj.chemical_shift{1,lo};%&&&&&&&&&&&&&&&&&&
%                         simple_data.main_cs_table_lin(  1,filled_by_other_compounds_cs+counter_spins)=assignment_obj.chemical_shift{1,lo};%&&&&&&&&&&&&&&&&&&
%                         simple_data.main_cs_table(  looop_over_compounds,counter_spins)=assignment_obj.chemical_shift{1,lo};%
%                         simple_data.lookup_cs_table2(looop_over_compounds,0*filled_by_other_compounds_cs+counter_spins)=lo;
%                         simple_data.lookup_cs_table_lin2                 (1,filled_by_other_compounds_cs+counter_spins)=lo;
%                          simple_data.lookup_cs_table(looop_over_compounds,0*filled_by_other_compounds_cs+counter_spins)=counter_spins;
%                          simple_data.lookup_cs_table_lin                 (1,filled_by_other_compounds_cs+counter_spins)=counter_spins;
%                         simple_data.lookup_compound_table2               (1,filled_by_other_compounds_cs+counter_spins)=looop_over_compounds;
%                         simple_data.label{                                 filled_by_other_compounds_cs+counter_spins}=label;
                        counter_spins=counter_spins+1;
                        
                    end
                    counter_param_spin=counter_param_spin+1;
                    if keep_label_core(1,counter_all_lables_including_skiped_ones)
                        counter_param_spin_returned=counter_param_spin_returned+1;
                    end
                    %add entry in list of symetrical spins for spinache
                    if how_many_atoms>1
                        bas.sym_group{1,counter_number_sym_group_magn_equ_spin}=['S' num2str(how_many_atoms)];
                        bas.sym_spins{1,counter_number_sym_group_magn_equ_spin}=loc_list_atom_number;
                        disp(['Introduced a sym group "S' num2str(how_many_atoms) '" in Spinach structure bas for spin labeles:<' labcur{1,1} '>.' ])
                        
                        counter_number_sym_group_magn_equ_spin=counter_number_sym_group_magn_equ_spin+1;
                        % bas.sym_spins={[14 15 16],[17 18 19],[20 21 22]};
                    end
                else
                    simple_data.main_cs_table(looop_over_compounds,lo)=0;
                 %%%  simple_data.main_cs_table2(looop_over_compounds,counter_spins)=0;
                end
                toto_check(1,counter_all_lables_including_skiped_ones)=1;
            else
                toto_check(1,counter_all_lables_including_skiped_ones)=0;
            end
            
            counter_all_lables_including_skiped_ones=counter_all_lables_including_skiped_ones+1;
            
        end
        filled_large=filled_large+size(assignment_obj.atom_num,1);
        
        if verbose
            if exist('list_diff_chem_shift','var')
                for lo=1:size(list_diff_chem_shift,1)
                    for l2=1:size(list_diff_chem_shift,2)
                        tm=list_diff_chem_shift(lo,l2);
                        if tm~=0
                            tmsc=assignment_obj.label_signal(1,list_label_number(lo));
                            tms=tmsc{1,1};
                            if verbose>2
                                disp(['Chemical shift number ' num2str(filled_by_other_compounds_chemical_shifts+tm) ' involves label number ' num2str(lo) ':' tms ' compound: ' num2str(looop_over_compounds)])
                            end
                        end
                    end
                end
            end
        end
        % scalar couplings
        %initialize
        a1=size(inter.coupling.scalar,1);
        a2=size(inter.coupling.scalar,2);
        for ll1=(a1+1):(a1+counter_spins-1)
            for ll2=(a2+1):(a2+counter_spins-1)
                inter.coupling.scalar{ll1,ll2}=[];
            end
        end
        for lop=1:size(j_obj.J,2)% loop over couplings found in J tag
            J=j_obj.J{1,lop};
            if abs(J)>min_j
                cur_atom_lab=j_obj.label_signal1{1,lop};
                n1=get_list_of_spin_number(assignment_obj,cur_atom_lab);%convert label into number
                list1=[];
                if keep_label(1,n1) % in the list of lables to keep
                    for lo=1:size(list_diff_chem_shift,1)
                        for l2=1:size(list_diff_chem_shift,2)
                            tm=list_diff_chem_shift(lo,l2);
                            if tm~=0
                                k=list_label_number(lo);
                                if k==n1
                                    % li1=[li1 lo];
                                    list1=[list1 tm];
                                    tmsc=assignment_obj.label_signal(1,k);
                                    tms=tmsc{1,1};
                                    %  disp(['First partner of couplings: ' num2str(tm) ' involves label number ' num2str(lo) ':' tms])
                                end
                            end
                        end
                    end
                end
                cur_atom_lab=j_obj.label_signal2{1,lop};
                n2=get_list_of_spin_number(assignment_obj,cur_atom_lab);
                list2=[];
                if keep_label(1,n2) % in the list of lables to keep
                    for lo=1:size(list_diff_chem_shift,1)
                        for l2=1:size(list_diff_chem_shift,2)
                            tm=list_diff_chem_shift(lo,l2);
                            if tm~=0
                                k=list_label_number(lo);
                                if k==n2
                                    list2=[list2 tm];
                                    tmsc=assignment_obj.label_signal(1,k);
                                    tms=tmsc{1,1};
                                    %  disp(['Second partner of couplings: ' num2str(tm) ' involves label number ' num2str(lo) ':' tms])
                                end
                            end
                        end
                    end
                end
                
                %  simple_data.main_j_table(1,filled_by_other_compounds_j+lop)=J;%keep_label_core
                simple_data.main_j_table(looop_over_compounds,lop)=J;%keep_label_core
                
                OK1 =  keep_label(1,n1) && keep_label(1,n2) ;% they are both in the list...
                OK2 = (keep_label(1,n1) && keep_label_core(1,n2) ) || ( keep_label_core(1,n1) && keep_label(1,n2) ) ;% they are both in the list and one of them in the core list
                %         if OK1 && (~OK2)
                %             disp('coupling outside the core')
                %         end
                if ( OK2 ) % they are both in the list and one of them in the core list
                    
                    
                    
                    %                          simple_data.main_j_table2d(ilookup(1,n1),ilookup(1,n2))=J;
                    %                          simple_data.main_j_table2d(ilookup(1,n2),ilookup(1,n1))=J;
                    
                    for ll1=1:size(list1,2)
                        list_diff_J1(filled_by_other_compounds_j+counter_param_coup,ll1)=list1(1,ll1);
                    end
                    for ll2=1:size(list2,2)
                        list_diff_J2(filled_by_other_compounds_j+counter_param_coup,ll2)=list2(1,ll2);
                    end
                    lookup_j_atom(filled_by_other_compounds_j+counter_param_coup)=counter_param_coup;
                    lookup_j_compound(filled_by_other_compounds_j+counter_param_coup)=looop_over_compounds;
                    
                    
                    
                    
                    if J~=0
                        for ll1=list1
                            for ll2=list2
                                if ll1<ll2
                                    a=ll1;
                                    b=ll2;
                                else
                                    a=ll2;
                                    b=ll1;
                                end
                                inter.coupling.scalar{a+filled_by_other_compounds_cs2,b+filled_by_other_compounds_cs2}=J;
                                if J>0
                                    diff_cs=(inter.zeeman.scalar{1,a} - inter.zeeman.scalar{1,b});
                                    if abs((larmor*diff_cs)/J)>10 % 10
                                        inter.coupling.strength{filled_by_other_compounds_cs2+a,filled_by_other_compounds_cs2+b}='zz';
                                        inter.coupling.strength{filled_by_other_compounds_cs2+b,filled_by_other_compounds_cs2+a}='zz';
                                    end
                                end
                            end
                        end
                    end
                    if J~=0
                        for ll1=list1
                            for ll2=list2
                                simple_data.main_j_table2d(looop_over_compounds,0+ll1,0+ll2)=J;
                                simple_data.main_j_table2d(looop_over_compounds,0+ll2,0+ll1)=J;
                            end
                        end
                    end
                    
                    simple_data.lookup_J_table(looop_over_compounds,0*filled_by_other_compounds_j+counter_param_coup)=lop;
                    %   simple_data.lookup_J_table_compound(1,filled_by_other_compounds_j+counter_param_coup)=looop_over_compounds;
                    counter_param_coup=counter_param_coup+1;
                end
            end
        end
        if verbose
            if exist('list_diff_J1','var')
                
                for lo=1:size(list_diff_J1,1)
                    for l2=1:size(list_diff_J1,2)
                        tm=list_diff_J1(lo,l2);
                        if tm~=0
                            %                 tmsc=assignment_obj.label_signal(1,list_label_number(tm));
                            %                 tms=tmsc{1,1};
                            if verbose>2
                                disp(['Coupling number ' num2str(lo) ' involves chemical shift number ' num2str(tm) ])
                            end
                        end
                    end
                    for l2=1:size(list_diff_J2,2)
                        tm=list_diff_J2(lo,l2);
                        if tm~=0
                            %                 tmsc=assignment_obj.label_signal(1,list_label_number(tm));
                            %                 tms=tmsc{1,1};
                            if verbose>2
                                
                                disp(['with chemical shift number ' num2str(tm) ])
                            end
                        end
                    end
                    if verbose>2
                        disp(['----------------------'])
                    end
                    
                end
            end
        end
        %tot.dfs=struct()
        
        if verbose
            inter.zeeman.scalar
            inter.coupling.scalar
            list_diff_chem_shift
            list_diff_J1
            list_diff_J2
        end
        
        % % % make list of H, separately (one line per H).
        % % % make list wht index increasing only when different labels (or forcae
        % % % equivalent ...) for the fitting to adjust only the desired ones ...
        % % 1.3, 1
        % % 1.3, 1
        % % 1.3, 1
        % % 2.5, 2
        % % 2.5, 2
        % % 4.8, 3
        % % ponters:
        % % 1: 1,2,3
        % % 2: 4,5
        % % 3: 6
        % % inter.zeeman.scalar={h_shift_1 c_shift_1 c_shift_2 h_shift_2};
        %For J...
        
        
        
        
        
        % % Absorb parameters
        % c_shift_1=params(1);
        % c_shift_2=params(2);
        % h_shift_1=params(3);
        % h_shift_2=params(4);
        % j_cc=params(5);
        % j_ch_far=params(6);
        % j_ch_near=params(7);
        % j_hh=params(8);
        % a_h=params(9);
        % a_c=params(10);
        % lw_h=params(11);
        % lw_c=params(12);
        %
        % % Spin system
        % sys.isotopes={'1H','13C','13C','1H'};
        %
        % % Magnet field
        % sys.magnet=2*pi*500.101412e6/spin('1H');
        %
        % % Chemical shifts
        % inter.zeeman.scalar={h_shift_1 c_shift_1 c_shift_2 h_shift_2};
        %
        % % Scalar couplings
        % inter.coupling.scalar={0.0   j_ch_near   j_ch_far  j_hh;
        %                        0.0   0.0         j_cc      j_ch_far;
        %                        0.0   0.0         0.0       j_ch_near;
        %                        0.0   0.0         0.0       0.0};
        %
        %     simple_data.start_pos_compound_cs(looop_over_compounds)=filled_by_other_compounds_chemical_shifts;
        %     simple_data.start_pos_compound_j(looop_over_compounds)=filled_by_other_compounds_j;
        filled_by_other_compounds_chemical_shifts=size(list_diff_chem_shift,1);
        filled_by_other_compounds_j=size(list_diff_J1,1);
        filled_by_other_compounds_cs2=filled_by_other_compounds_chemical_shifts;
        counter_spins
        filled_by_other_compounds_j=filled_by_other_compounds_cs;%%%TTTT
        filled_by_other_compounds_cs=counter_spins-1;
    end
    spinach_system.sys=sys;
    spinach_system.inter=inter;
    spinach_system.bas=bas;
    
    
    if verbose==3
        disp(['Compound : ' num2str(looop_over_compounds) ' out of ' num2str(size(item_num,2)) '========================================' ])
    end
end
optim_param.lookup_cs_atom=lookup_cs_atom;
optim_param.lookup_cs_compound=lookup_cs_compound;
optim_param.lookup_j_atom=lookup_j_atom;
optim_param.lookup_j_compound=lookup_j_compound;
optim_param.list_diff_J1=list_diff_J1;
optim_param.list_diff_J2=list_diff_J2;
optim_param.list_diff_chem_shift=list_diff_chem_shift_returned;
end
function                 n=get_list_of_spin_number(assignment_obj,cur_atom_lab)
%remove space in case start with one...
if  cur_atom_lab(1,1)==32
    cur_atom_lab=cur_atom_lab(1,2:end);
    % warning('Remove a space in label... in J_tag')
end
n=0;
for l=1:size(assignment_obj.label_signal,2)
    lab=assignment_obj.label_signal{1,l};
    %   disp(['test ' cur_atom_lab '  label...' lab])
    
    if strcmp(cur_atom_lab,lab)
        n=l;
        break
    end
end
if n==0
    for l=1:size(assignment_obj.label_signal,2)
        lab=assignment_obj.label_signal{1,l};
        disp(['<' lab '>'])
    end
    warning(['Could not find <' cur_atom_lab '>  label... in the list above... check the NMREDATA_ASSIGNMENT or NMREDATA_J tags'])
end
end