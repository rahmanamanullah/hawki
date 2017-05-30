Optimized near infrared imaging sky subtraction.

Sky images are built based on consecutive near infrared dithered images, where
the tidal sky variation between images is expected to be small.  A running 
mean/median over time is used for estimating the sky value in each pixel, and
one sky frame and mask is saved for each image.  

SExtractor (https://www.astromatic.net/software/sextractor) can be used to create
object masks for each exposure which is then taken into account when building the
sky frame.

This project was initiated with repeated observations of the same field in mind,
a strategy that is typical for transient searches.  A deep object mask can be
be continuously updated when the same field is observed repeatedly.  Previous
observations are then rereduced using the improved object mask.  The result is
a much cleaner image that allows searches for low signal-to-noise transients.

This code is using CCFits (https://heasarc.gsfc.nasa.gov/fitsio/ccfits/) to
handle the FITS IO.

The code was used for the analyses presented in:

- Amanullah et al. (2011) :    https://arxiv.org/abs/1109.4740
- Petrushevska et al. (2016) : https://arxiv.org/abs/1607.01617
- Petrushevska et al. (2017) :
