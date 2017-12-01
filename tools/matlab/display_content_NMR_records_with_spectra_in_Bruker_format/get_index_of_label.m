function  [ind]=get_index_of_label(object,lab)
ind=0;
%first run for explicit atoms
for lo_over_labels=1:size(object.label_signal,2)
    if strcmp(object.label_signal{1,lo_over_labels},lab)
        ind=lo_over_labels;
        break
    end
end
   