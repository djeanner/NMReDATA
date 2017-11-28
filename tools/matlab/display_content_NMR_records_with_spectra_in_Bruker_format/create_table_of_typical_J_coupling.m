function [J_structure,typ_coupling]=create_table_of_typical_J_coupling(min_value_to_includeJ)
typ_coupling=zeros(92,92,4);%for 92 pairs of elements / 1J-4J
typ_coupling(1,1,1)=480;%1J HH
typ_coupling(1,1,2)=16;%1J HH
typ_coupling(1,1,3)=3;%1J HH
typ_coupling(1,1,4)=0.9;%1J HH
typ_coupling(1,6,1)=160;%1J CH
typ_coupling(1,6,2)=5;%2J CH
typ_coupling(1,6,3)=3;%3J CH
%add to 31P, 19F...

%1H
J_structure{1}=[1/2 1/2];% coupling to 1h for element 1, only spin 1/2 taken into account, doublet.

%13C
if min_value_to_includeJ>1
    J_structure{6}=[1];% When neglecting 1% 13C
else
    J_amplitudes{6}=[1/2 99 1/2]/100;% coupling to 13C: main signal 99%, 1% 13C satellites, spin 1/2,
end

% Sn
J_structure{50}=[100]/100;% for Sn 116,118,120,122,124 S0(117Sn 7.68% 119Sn 8.59% ... gamma ?
if min_value_to_includeJ<8.59
    J_structure{50}=[ 8.59/2 (14.54+24.22+32.58+4.63+5.79)  8.59/2 ]/100;% for Sn 116,118,120,122,124 S0(117Sn 7.68% 119Sn 8.59% ... gamma ?
end
if min_value_to_includeJ< 7.68
    J_structure{50}=[7.68/2 8.59/2 (14.54+24.22+32.58+4.63+5.79)  8.59/2 7.68/2]/100;% for Sn 116,118,120,122,124 S0(117Sn 7.68% 119Sn 8.59% ... gamma ?
end

%add Si

end