clear all

%file version 1.1 (that is with backslash)
nmredata=check_nmredata_sdf_file('./NMR_Records_folder/HAP_benzo(a)pyrene_assignments/compound1.nmredata_V1.1.sdf');

%file version 1.0 (that is with backslash)
%nmredata=check_nmredata_sdf_file('./NMR_Records_folder/HAP_benzo(a)pyrene_assignments/compound1.nmredata.sdf');

export_to_acd(nmredata);