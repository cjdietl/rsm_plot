# rsm_plot

Converting Angular Space of a 1D-Detector to a Reciprocal Space Map (RSM)

During the analysis of epitaxial thin films, one of the first things to study is the strain relative to the substrate. This is best  done by visualizing the peak positions of the thin film in relation to the peaks of the substrate the film was grown on. This is best done by a so called 'Reciprocal Space Map' or RSM for short. Here, one selects in-plane Bragg reflections of both substrate and a thin film and rasters them with an X-ray diffractometer. In this way, the amount of strain can be quantified. See e.g. one of my publications ( http://aip.scitation.org/doi/10.1063/1.5007680) for some examples.

This particular MatLab-script takes the pixels of a strip-detector (Mython) and combines them with a so SPEC-File. A SPEC-File is a data file generated by SPEC, a very common programm used by many X-ray spectroscopists to take data and to navigate reciprocal space. After applying some linear algebra, the angular information from the file is converted to reciprocal space. This is done in rsm_map_fourc.m . This generated a HDF-file which containts all reciprocal space coordinates of all pixels. However, these coordinates are irregularly spaced and have to be interpolated and gridded in order to make an image. This happens in create_gridded_rsm_view.m .






