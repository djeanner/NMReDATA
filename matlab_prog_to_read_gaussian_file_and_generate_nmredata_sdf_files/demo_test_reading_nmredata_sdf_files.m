% demo_test_reading_nmredata_sdf_files.m  check .nmredata.sdf files
%% select .nmredata.sdf file to check
file_name='./nmredata_sdf_demo_files/benzoapyrene.sdf';
file_name='./nmredata_sdf_demo_files/benzocchrysene.sdf';
file_name='./nmredata_sdf_demo_files/androsten.sdf';

%% call main reading function
% the main function reading and testing nmredata sdf files is:
verbosity=1;%level of verbosity
[super_obj, returned_value, text_of_the_problem]=check_nmredata_sdf_file(file_name,verbosity);

%% display results
if returned_value
    disp('Object OK!')
    super_obj
else
        disp('Reading .nmredata.sdf files failed!')
        disp(text_of_the_problem)
end
    
    