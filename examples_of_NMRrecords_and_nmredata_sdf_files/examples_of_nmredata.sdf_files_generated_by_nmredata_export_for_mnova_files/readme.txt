A set of nmredata.sdf and NMR records were generated from nova files using a Mnova script. They are just for illustration and test purpose. They were generated using a export script written by a Damien Jeannerat. Dont hesitate to contact me (damien.jeannerat@unige.ch) if you have any question... This dropbox folder is related to the github repository

Link to the script: https://github.com/djeanner/NMReDATA/tree/master/tools/Generate_NMR_record_from_Mnova_files/Mnova_script_to_generate_of_nmredata_file

Related source of data:

1) Only the nmredata.sdf files can be found in GitHub because of space limitation reasons.
https://github.com/djeanner/NMReDATA/tree/master/examples_of_NMRrecords_and_nmredata_sdf_files

2) full NMR records and corresponding Mnova files can be found in this dropbox:
https://www.dropbox.com/sh/ma8v25g15wylfj4/AAA4xWi5w9yQv5RBLr6oDHila?dl=0

Comments about the samples:
a-b-glucose contains more than one assigned compound
3,5.3,5-Bis(trifluoromethyl)aniline includes multiplet in the 13C spectrum because of coupling with 19F. It includes HOESY spectra. Some versions include the assignment of the HOESY spectra.
Caryophyllene_oxide_full_assignments includes two non-implicit hydrogens.
Menthol contains a lot of non-equivalent H on CH2 and non-classical spectra.

Important note: Only “true” assignment are extracted. When peaks from 2D spectra are linked to atoms of molecules, they are not (yet) taken into account. (for example in Ice_tea_lemon_partial_assignments_HSQC_assignement_of_peaks.mnova)