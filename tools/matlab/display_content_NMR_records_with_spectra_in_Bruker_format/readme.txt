This matlab/octave program is under developments.
As Nov. 9th, this is what it does:
- Reads the nmredata.sdf file in unzip NMR record
- Reads the nmr spectra (in Bruker format) of for the data in the nmredata.sdf file
- It plots figures of the 1D and 2D spectra (in blue) and include the NMReDATA on top (in black). 
- It generates .pdf of the figures.


When octave (a free clone of matlab) is installed, the demo can be run as follows:
octave  --no-gui --no-window-system demo_display_NMR_record.m 
