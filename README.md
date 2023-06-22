# MESOX

The keywords in the headers of MAROON-X data is not compatible with the required format to run ESO's molecfit on it (https://www.eso.org/sci/software/pipelines/skytools/molecfit). 
This tool serves as a translator to enable smooth application when running tayph. You can find the tayph-repo here: https://github.com/Hoeijmakers/tayph.

Note that this code requires that you have access to all the raw files of from the MX observations. Reach out to your lovely night astronomer for it!
This code allows you to use molecfit within the framework of tayph. 

If you don't want to use molecfit with tayph, you can still use the .par file and the inclusion regions. You can also stitch the spectra (in molecfit_on_MX.py), but molecfit outside of tayph requires a different format.
