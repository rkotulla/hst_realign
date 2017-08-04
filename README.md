# HST Re-align

This repository is a small tool to fix pointing offsets in Hubble Space 
Telescope (HST) data. However, the tool is also applicable to any other FITS
data as long as the World Coordinate System (WCS) is described using the usual 
CRVAL/CRPIX/CD header keywords.
 
## How it works

Upon starting the program automatically opens the specified (reference-) image 
 in ds9, and overlays a 2MASS point source catalog. The user is then asked to 
 select pairs of coordinates by clicking th mouse in the ds9 window. Each pair 
 starts with the position of a well-detected source in the image, followed by 
 its counterpart from the 2MASS catalog. In other words, figure out what star 
 corresponds to what source, then first click on the star, then on the green 
 circle marking the position in 2MASS.
  
Once this is repeated for a small number of stars (at least three, spaced out 
across the image for optimal precision), the code then computes a rough 
estimate of the corrected WCS, and displays the best-guess image in ds9, 
again overlayed with the 2MASS source catalog.

In the next step, the code computes a full source catalog using Astromatic's 
Sextractor to extract accurate position for all sources in the field. It also 
downloads a list of calibrated sources from the GAIA survey. Both catalogs are
then matched (using a matching radius of 1.0"), providing a new and improved 
set of X/Y positions and corresponding sky-positions in Ra/Dec that allow to
optimize the WCS system. This step of matching catalogs and optimizing the WCS
is then iterated to provide the best possible solution using the largest 
possible number of sources.

in a final step, the new WCS system is then applied to all input frames, with 
the corrected files being written back to disk using a suffix (default is 
wcsfix) inserted in front of the .fits extension.



## How to run it