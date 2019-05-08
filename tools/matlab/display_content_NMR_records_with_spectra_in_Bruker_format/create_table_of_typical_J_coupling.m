function [J_structure,typ_coupling]=create_table_of_typical_J_coupling(min_abundance_take_into_account,min_value_to_includeJ)
if nargin==0
    min_abundance_take_into_account=0.5/100;% below 1/100 will take into account 13C sat unit
    min_abundance_take_into_account=0.01/100;% Deuterium sat..
        min_abundance_take_into_account=1.5/100;% standard...
    min_value_to_includeJ=0;

end

verbose=1;
if verbose
   disp( 'Construction of data for multiplet reconstruction (structure and values)')
end
istope_prop(1,1,:)=[100  1/2 99.9845/100];%1H
istope_prop(1,2,:)=[15.350609  1 0.0155/100];%2H (D)
%istope_prop(1,3,:)=[106.663974 1/2 3e-16/00];%3H
istope_prop(6,13,:)=[25.145020 1/2 1.108/100];%1H

istope_prop(50,117,:)=[35.632259  1/2 7.68/100];%117Sn
istope_prop(50,119,:)=[32.718749 1/2 8.59/100];%119Sn
istope_prop(9,19,:)=[94.094011 1/2 100/100];%19F
istope_prop(15,31,:)=[40.480742 1/2 100/100];%31P
istope_prop(5,10,:)=[10.743658 3 19.9/100];%10B
istope_prop(5,11,:)=[32.08397 3/2 80.1/100];%11B

%table of typical coupling
typ_coupling=zeros(92,92,4);%for 92 pairs of elements / 1J-4J
typ_coupling(1,1,1)=480;%1J HH ?? Not sure should be useless
typ_coupling(1,1,2)=16;%2J HCH
typ_coupling(1,1,3)=3;%3J HCCH
typ_coupling(1,1,4)=0.9;%1J HH
typ_coupling(1,6,1)=160;%1J CH
typ_coupling(1,6,2)=5;%2J CCH
typ_coupling(1,6,3)=3;%3J CCCH
typ_coupling(1,6,4)=0.5;%4J CCCCH
typ_coupling(6,50,1)=324.3782;%1J CSn (12.2741-10.1259)*151=324.3782 %119Sn
typ_coupling(6,50,2)=1.9560;%1J CSn (29.1743-29.0439)*151=1.9560 %119Sn
typ_coupling(1,6,1)=65;%1J HB
typ_coupling(6,12,1)=75;%1J BC %from DOI: 10.1021/ic50145a028
typ_coupling(1,9,2)=50;%2J HF 
typ_coupling(1,9,3)=13;%3J HF 
typ_coupling(6,9,1)=272;%1J CF 
typ_coupling(6,9,2)=33;%2J CF 
typ_coupling(6,9,3)=4;%3J CF 
typ_coupling=typ_coupling.*(typ_coupling>min_value_to_includeJ);
%% add to 31P, 19F...


% This is obsolete
% %1H
% J_structure{1}=[1/2 1/2];% coupling to 1h for element 1, only spin 1/2 taken into account, doublet.
% 
% %13C
% if min_value_to_includeJ>1
%     J_structure{6}=[1];% When neglecting 1% 13C
% else
%     J_structure{6}=[1/2 99 1/2]/100;% coupling to 13C: main signal 99%, 1% 13C satellites, spin 1/2,
% end
% 
% % Sn
% J_structure{50}=[100]/100;% for Sn 116,118,120,122,124 S0(117Sn 7.68% 119Sn 8.59% ... gamma ?
% if min_value_to_includeJ<8.59
%     J_structure{50}=[ 8.59/2 (14.54+24.22+32.58+4.63+5.79)  8.59/2 ]/100;% for Sn 116,118,120,122,124 S0(117Sn(35.632259) 7.68% 119Sn(32.718749) 8.59% ... gamma ?
% end
% if min_value_to_includeJ< 7.68
%     J_structure{50}=[7.68/2 8.59/2 (14.54+24.22+32.58+4.63+5.79)  8.59/2 7.68/2]/100;% for Sn 116,118,120,122,124 S0(117Sn 7.68% 119Sn 8.59% ... gamma ?
% end
% 
% %add Si
% J_structure2=J_structure;


for loop_el=1:size(istope_prop,1)
    if  sum(istope_prop(loop_el,:,2))>0
        J_structure{loop_el}=[];
        J_amplitudes{loop_el}=[];
        others_abundance=1;
        % most abundant isotope
        tmp_for_sort=istope_prop(loop_el,:,3);
        [~,tmp_for_sort]=sort(tmp_for_sort,'descend');
        ref_gam=istope_prop(loop_el,tmp_for_sort(1,1),1);
        for loop_is=1:size( istope_prop,2)
            if istope_prop(loop_el,loop_is,3)>min_abundance_take_into_account
                if istope_prop(loop_el,loop_is,2)>0
                    nb_line=(istope_prop(loop_el,loop_is,2)*2)+1;
                    scaling=istope_prop(loop_el,loop_is,1)/ref_gam;
                    for loop_over_transition=1:nb_line
                        J_structure{loop_el}=[J_structure{loop_el} istope_prop(loop_el,loop_is,3)/nb_line];
                        J_amplitudes{loop_el}=[ J_amplitudes{loop_el} scaling];
                        others_abundance=others_abundance-istope_prop(loop_el,loop_is,3)/nb_line;
                    end
                end
            end
        end
        if others_abundance>min_abundance_take_into_account
            J_structure{loop_el}=[J_structure{loop_el} others_abundance];
            J_amplitudes{loop_el}=[ J_amplitudes{loop_el} 0];
        end
        %normalize
        J_structure{loop_el}=J_structure{loop_el}/(sum(J_structure{loop_el}));
    end
end

if verbose
for loop_el=1:size(istope_prop,1)
    tmp=J_structure{loop_el};
    if size(tmp,1)>0
        disp(['Element:' num2str(loop_el) ' structure/amplitudes:' num2str(tmp)])
    end
    tmp=J_amplitudes{loop_el};
    if size(tmp,1)>0
        disp(['Element:' num2str(loop_el) ' rel.split:' num2str(tmp)])
    end
end
end


end