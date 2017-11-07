function [ ret_nois log_ref]=determine_noise_level(data)

%if type_of_selection==1

spectrum=data.spectrum;
%function [topg ret_nois cont_level_list]=plot_noise_of_spectrum_2(base_path,experiment_name,expno,procno,no_clear,sim,vert_scale,how_much_higher,factor_between_level_contours,how_much_higher2)
%function [topg ret_nois cont_level_list]=plot_noise_of_spectrum_2(base_path,experiment_name,expno,procno,no_clear,sim,vert_scale,how_much_higher,factor_between_level_contours,how_much_higher2)

ret_nois=0;


log_ref=         (spectrum(2:end-1,2:end-1)>spectrum(3:end  ,2:end-1));
log_ref=log_ref.*(spectrum(2:end-1,2:end-1)>spectrum(1:end-2,2:end-1));
log_ref=log_ref.*(spectrum(2:end-1,2:end-1)>spectrum(2:end-1,1:end-2));
log_ref=log_ref.*(spectrum(2:end-1,2:end-1)>spectrum(2:end-1,3:end  ));

log_ref=log_ref.*(spectrum(2:end-1,2:end-1)>spectrum(3:end-0,3:end-0));
log_ref=log_ref.*(spectrum(2:end-1,2:end-1)>spectrum(1:end-2,3:end-0));
log_ref=log_ref.*(spectrum(2:end-1,2:end-1)>spectrum(3:end-0,1:end-2));
log_ref=log_ref.*(spectrum(2:end-1,2:end-1)>spectrum(1:end-2,1:end-2));

log_ref=log_ref.*(spectrum(2:end-1,2:end-1)>0);
log_ref=log_ref.* spectrum(2:end-1,2:end-1);
log_ref=reshape(log_ref,size(log_ref,1)*size(log_ref,2),1);
pea_list_find=find(log_ref~=0);
log_ref=log_ref(pea_list_find,1);
log_ref=sort(log_ref,'descend');
if size(log_ref,1)>0,
    ret_nois=log_ref(round(size(log_ref,1)/2),1);
else
    ret_nois=0;
    
    
end
end