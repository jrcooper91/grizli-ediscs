
# coding: utf-8

# The `Grizli` pipeline allows you to fully reduce a given set of HST grism observations with essentially two steps:
# 
# * Run an archive query with [`hsaquery`](https://github.com/gbrammer/esa-hsaquery)
# 
# * Process the associations found with the query with `grizli.auto_script.go`.
# 
# Here, "association" usually simply means "any Hubble exposures that overlap" and doesn't require that all observations were taken with the same observing program, instrument, grism position angle, epoch, filter, etc.  The code does all of the exposure-level book-keeping and the products are drizzled image mosaics, extracted 1D and 2D grism spectra and fits to the spectra.
# 
# **NB**: The pipeline works fine with just imaging and no grism exposures!

# In[7]:


#cell 1
#get_ipython().run_line_magic('matplotlib', 'inline')


# In[8]:


#cell 2
import glob
import time
import os

import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Image

import astropy.io.fits as pyfits
import drizzlepac

import grizli
from grizli.pipeline import auto_script 
from grizli import utils
from grizli import fitting
from grizli import multifit #original line, replaced by cell above 

utils.set_warnings()
print('\n Grizli version: ', grizli.__version__)


# In[9]:


#cell 3
#os.chdir('/Users/brammer/3DHST/Spectra/Work/Grizli/Demo-18.05.22/')
HOME_PATH = os.getcwd()
print('HOME_PATH = ', HOME_PATH)


# ## Query the HST archive ##
# 
# The `hsaquery` module can be used to programaticaly query the HST archive and find exposures from different programs (and instruments) that overlap on the sky.  The example below is tailored for a single pointing from a single program, but the query parameters can be expanded to search much more broadly for archival data.

# In[10]:


#cell 4
### Generate a query for the WFC3/ERS grism data

## !! new query tools since ESA database changed in summer 2018
# https://github.com/gbrammer/mastquery
from mastquery import query, overlaps

# "parent" query is grism exposures in GO-11359.  Can also query the archive on position with
# box=[ra, dec, radius_in_arcmin]
parent = query.run_query(box=None, proposal_id=[12945], instruments=['WFC3/IR', 'ACS/WFC'], 
                         filters=['G102','G141'])

# ### "overlap" query finds anything that overlaps with the exposures 
# ### in the parent query
# extra = query.DEFAULT_EXTRA # ignore calibrations, etc.

# ## To match *just* the grism visits, add, e.g., the following:
# extra += ["TARGET.TARGET_NAME LIKE 'WFC3-ERSII-G01'"]

tabs = overlaps.find_overlaps(parent, buffer_arcmin=0.01, 
                              filters=['F105W', 'F814W','G102'], 
                              proposal_id=[12945], instruments=['WFC3/IR','WFC3/UVIS','ACS/WFC']) 
                              #,extra={'target_name':'CL1059-12.0'}, close=False)


# In[11]:


#cell 5
# Summary of the tables you just generated
foot_files = glob.glob('j[02]*footprint.fits')
print('Footprint files: ', foot_files)

print('\n# id            ra         dec        e(b-v)   filters')
for tab in tabs:
    print('{0}  {1:.5f}  {2:.5f}   {3:.4f}   {4}'.format(tab.meta['NAME'], tab.meta['RA'], 
                                                 tab.meta['DEC'], tab.meta['MW_EBV'],
                                                  ','.join(np.unique(tab['filter']))))


# In[12]:


#cell 6
#os.chdir('/Users/brammer/3DHST/Spectra/Work/Grizli/Demo-18.05.22/')
HOME_PATH = os.getcwd()
print('HOME_PATH = ', HOME_PATH)


# # - Pipeline processing - #
# 
# ** In principle, all of the steps outlined below can be executed with a single call to** `auto_script.go`, from fetching the data to extracting spectra and performing the redshift / line fits.  The processing steps been broken out individually here to show the processing at each step.
# 
# ** The same pipeline can be used to process imaging-only fields.**  Simply run the queries as above to find the imaging exposures you want to processes and run everything the same way.  The pipeline steps related to the grism exposures will simply be skipped.

# In[14]:


#cell 7
# Do everything for the query from fetching the data to generating the contamination model
HOME_PATH = os.getcwd()
print('HOME_PATH = ', HOME_PATH)

"CHANGE FILE DIRECTORY HERE"
root = 'j122816m1132'
IS_PARALLEL = False # Set to True for parallel programs like WISPS

if False:
    # This line would do everything below
    auto_script.go(root=root, maglim=[19,21], HOME_PATH=HOME_PATH, reprocess_parallel=True, 
                   s3_sync='cp', gaia_by_date=True, is_parallel_field=IS_PARALLEL, 
                   run_fit=False, only_preprocess=True, run_extractions=False)


# # - Individual steps - #
# 
# ## Fetch data from the HST archive ##
# `Grizli` can automatically fetch HST data from the ESA Hubble Science archive (and, optionally, the Amazon S3 bucket).  The `fetch_files` script fetches the exposures listed in the archive query above.  It also fetches associated WFC3/IR persistence products from the persistence database.
# 
# The first time you run the script, a lot more information will be printed to the screen as the exposures are retrieved and the script runs the reprocessing code to flatten the IR backgrounds.  Below the "skip" message simply indicate that files have already been downloaded.

# In[15]:


#cell 8 THIS TAKES A LONG TIME
### Fetch data, reprocess WFC3/IR for backgrounds, fetch WFC3/IR persistence productss

# If s3_sync, then sync from the Hubble Amazon S3 bucket with awscli, 
# otherwise get from the ESA archive.
os.chdir(HOME_PATH)

import grizli.pipeline
from grizli.pipeline import auto_script
# Is awscli available and connected? 
s3_status = os.system('aws s3 ls s3://stpubdata --request-payer requester')
if s3_status == 0:
    s3_sync='cp'  # As of late October 2018, 's3 sync' not working with 'stpubdata'
else:
    s3_sync=False # Fetch from ESA archive
    
auto_script.fetch_files(field_root=root, HOME_PATH=HOME_PATH, remove_bad=True, 
                        reprocess_parallel=True, s3_sync=s3_sync)


# ## Parse visit associations ##
# `Grizli` builds its own associations based on anything it finds in the `RAW` directory.  Visits are usually defined in the exposure filenames.  For example, for the single exposure, `ib6o03ntq_flt.fits`, the characters `b6o` identify the observing program and the visit identifier  is `03`.  You can also build visits combining all exposures in a given filter taken at the same position angle, which can be useful for some programs executed in parallel where exposures taken at a similar time could have different visit IDs in the filename.  
# 
# **NB:** Generally one should process "visits" as groups of exposures in a given filter that were taken with a single guide star acquisition.  
# 
# The parsing script also associates grism exposures with corresponding direct images, based on the visit, exposure order and exposure footprints on the sky.

# In[16]:


#cell 11
# Demo combining by PA / filter.  

# Here it actually gets a bit confused because multiple F098M exposures 
# were taken at the same PA but shouldn't be associated with the grism exposures.
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))
visits, all_groups, info = auto_script.parse_visits(field_root=root, 
                                                    HOME_PATH=HOME_PATH, use_visit=True, 
                                                    combine_same_pa=True)

print('\n ====== \n')
for visit in visits:
    print('{0:30} {1:>2d}'.format(visit['product'], len(visit['files'])))


# In[17]:


#cell 12
######################
### Parse visit associations for most normal programs
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))
visits, all_groups, info = auto_script.parse_visits(field_root=root, 
                                                    HOME_PATH=HOME_PATH, use_visit=True, 
                                                    combine_same_pa=IS_PARALLEL)

print('\n ====== \n')
for visit in visits:
    print('{0:30} {1:>2d}'.format(visit['product'], len(visit['files'])))


# ## Master Pre-processing script: `grizli.prep.process_direct_grism_visit` ##
# 
# The `process_direct_grism_visit` script in [prep.py](https://github.com/gbrammer/grizli/blob/master/grizli/prep.py) provides one-stop-shopping for all of the preprocessing steps required.  This includes
# 
# * File handling (e.g., copying from `./RAW` to `./Prep/`)
# * Astrometric registration
# * Grism sky background subtraction & flat-fielding
# * Extract visit-level catalogs and segmentation images from the direct imaging
# 
# The products of the script for a given direct/grism pair are 
# 
# * Aligned, background-subtracted FLTs
# * Drizzled mosaics of direct & grism images
# 
# The script also runs on *imaging-only* visits, performing the background subtraction and astrometric alignment but skipping anything related to grism processing.
# 
# The `auto_script.preprocess` command below runs the processing script for the two direct/grism pairs of the ERS observations and for the overlapping imaging visits identified in the initial query.  It prints a bunch of information to the terminal, primarily from various runs of AstroDrizzle, and takes a few minutes to run per visit.  It only needs to be run once.
# 
# **NB** If you restart the pipeline after a previous run, it will skip preprocessing any visit where the file `{visit-product}_dr?_sci.fits` is found (i.e., the "Skip" messages below).  If you want to force reprocessing of a visit, delete that file.

# In[18]:


#####################
### Alignment & mosaics    
os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

# Alignment reference catalogs, searched in this order
catalogs = ['NSC', 'PS1','SDSS','GAIA','WISE'] 
# As of v0.8.0-8, can use the NOAO source catalog (NSC) here, which 
# is defined over much of the sky and appears to be well aligned to GAIA.  
# However, sometimes it's not clear how to apply the best quality control 
# to the NSC sources.  Here, for example, there seem to be a number of spurious 
# NSC sources that make the initial alignment RMS fairly high. 

# This script will do all the preprocessing of the grism *and* imaging visits 
# found in your archive query.
auto_script.preprocess(field_root=root, HOME_PATH=HOME_PATH, 
                       make_combined=False, catalogs=catalogs, use_visit=True)


# In[19]:


#get_ipython().system('ls wfc3*sci.fits # individual drizzled visits')


# In[20]:


# Results of the intra-visit alignment.  
# Should be small as these are just FGS drift on a single guide star
#get_ipython().system('ls *shifts.log')
print('')
#get_ipython().system('cat *shifts.log')


# In[21]:


# Show the alignment w.r.t the external NOAO Source Catalog (NSC)
#Image(filename = "./cl1138-11.0-c1b-05-122.0-f105w_wcs.png") 


# In[22]:


# Show the alignment of one HST visit to another, note tight 
# plot range compared to previous
#Image(filename = "./cl1138-11.0-c1b-05-122.0-g102_column.png") 


# In[23]:


# Check wcs.log files that a few objects were found, no large rotations
# and rms (second to last column) isn't too large

# Here, the F098M/G102 visit was the first processed and was aligned 
# to the NSC, with RMS~0.8 WFC3/IR pix.  Subsequent visits are aligned to 
# previously processed HST visits so that at least the relative HST astrometry
# is as good as possible.  Here, the F140W/G141 visit was run second and was 
# therefore aligned to F098M, resulting in much better precision than with the
# external catalog (RMS < 0.1 pix).

# Cat wcs.log files in order they were generated
#get_ipython().system('grep " 0 " `ls -ltr *wcs.log | awk \'{print $9}\'` | sed "s/  */ /g"')

# columns: 
# "visit"  0  xshift yshift rot scale rms N


# ### Alignment failures ###
# 
# The main failure mode of the `auto_script.preprocess` script is failure to compute a reliable alignment to the external reference.  This can happen, e.g., if there are not enough alignment sources (i.e., zero) within the field of view or if the original astrometry of the exposures obtained from the archive is offset from the reference by more than ~10 pixels.  This can almost always be remedied by running `grizli.pipeline.auto_script.manual_alignment` after the files have been fetched, which prompts the user to interactively mark sources in the image and reference catalog using DS9.

# In[24]:


if False: # Don't run
    catalogs = ['PS1','SDSS','GAIA','WISE']
    auto_script.manual_alignment(field_root=root, HOME_PATH=HOME_PATH, skip=True, 
                                 catalogs=catalogs, radius=15, visit_list=None)


# ### Grism sky subtraction ###
# 
# The grism sky backgrounds are subtracted using the "Master sky" images from [Brammer, Ryan, & Pirzkal 2015](http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2015-17.pdf) (available [here](http://www.stsci.edu/hst/wfc3/analysis/grism_obs/calibrations/wfc3_grism_master_sky.html)).  
# 
# `Grizli` ignores any wavelength dependence of the flat-field and applies a grey correction using the F140W (F105W) flat-field for the G141 (G102) grisms.
# 
# Residuals after subtracting the master sky images are typically of order 0.005 e-/s, just 0.5-1% overall background level.  They are removed by subtracting a column-average of the sky pixels in the grism exposures, and the processing script produces a diagnostic figure like the one shown below.  The source of the residuals is still unclear (e.g., perhaps spectra of objects near/below the detection limit).  Though they are usually well removed by the column average, they do make extracting continuum spectra of faint sources challenging.

# In[25]:


#from IPython.display import Image
#Image(filename = "./cl1138-11.0-c1b-05-122.0-g102_column.png", width=600, height=600)


# ### Fine alignment to GAIA DR2 ###
# 
# The initial visit alignment scripts often show small drifts such that the differen't visits don't perfectly overlap.  The script below performs an additional realignment to the visits internally and also to an external reference, usually GAIA DR2.

# In[26]:


# Fine alignment of the visits relative to each other and absolute to GAIA DR2
if len(glob.glob('{0}*fine.png'.format(root))) == 0:
    fine_catalogs = ['GAIA','PS1','SDSS','WISE']
    out = auto_script.fine_alignment(field_root=root, HOME_PATH=HOME_PATH, 
                                     min_overlap=0.2, stopme=False, ref_err=0.08, 
                                     catalogs=fine_catalogs, NITER=1, maglim=[17,23],
                                     shift_only=True, method='Powell', redrizzle=False, 
                                     radius=10, program_str=None, match_str=[], 
                                     gaia_by_date=True)

    # Update headers with the result from the fine alignment
    # Original FLTs are archived to FineBkup
    auto_script.update_wcs_headers_with_fine(root)
    
visits, res = np.load('{0}_fine.npy'.format(root))
shifts = res.x.reshape((-1,2))/10.
for i, visit in enumerate(visits):
    print('{0:35}  {1:6.2f}  {2:6.2f}'.format(visit['product'], shifts[i,0], shifts[i,1]))


# In[27]:


# Show the results of fine alignment.  
# Top panels are alignment between the visits.  + in the bottom panels are 
# residuals of the external reference, here GAIA DR2.
#
# Small drift between individual visits removed.  
# Fairly large GAIA offsets probably due to ~6 years between 
# WFC3/ERS and GAIA epoch 2015.5.
Image(filename='{0}_fine.png'.format(root))


# ## Make combined mosaics for each available filter ##
# 
# These are used to generate a photometric catalog and also for the direct image reference for the grism

# In[28]:


# Drizzle mosaics in each filter and combine all IR filters
combine_all_filters=True
if len(glob.glob('{0}-ir_dr?_sci.fits'.format(root))) == 0:

    ## Mosaic WCS
    wcs_ref_file = '{0}_wcs-ref.fits'.format(root)
    if not os.path.exists(wcs_ref_file):
        auto_script.make_reference_wcs(info, output=wcs_ref_file, 
                           filters=['G800L', 'G102', 'G141'], 
                           pad_reference=90, pixel_scale=None,
                           get_hdu=True)

    # All combined
    IR_filters = ['F105W', 'F110W', 'F125W', 'F140W', 'F160W', 
                  'F098M', 'F139M', 'F127M', 'F153M']

    optical_filters = ['F814W', 'F606W', 'F435W', 'F850LP', 'F702W', 'F555W', 'F438W', 'F475W', 'F625W', 'F775W', 'F225W', 'F275W', 'F300W', 'F390W']

    if combine_all_filters:
        auto_script.drizzle_overlaps(root, 
                                 filters=IR_filters+optical_filters, 
                                 min_nexp=1, 
                                 make_combined=True,
                                 ref_image=wcs_ref_file,
                                 drizzle_filters=False) 

    ## IR filters
    auto_script.drizzle_overlaps(root, filters=IR_filters, 
                                 min_nexp=1, 
                                 make_combined=(not combine_all_filters),
                                 ref_image=wcs_ref_file) 

    # Fill IR filter mosaics with scaled combined data so they can be used 
    # as grism reference
    auto_script.fill_filter_mosaics(root)

    ## Optical filters

    mosaics = glob.glob('{0}-ir_dr?_sci.fits'.format(root))

    auto_script.drizzle_overlaps(root, filters=optical_filters,
        make_combined=(len(mosaics) == 0), ref_image=wcs_ref_file,
        min_nexp=2) 


# In[29]:


#get_ipython().system('ls -1 j*_dr?_sci.fits')


# ## Generate a photometric catalog ##
# 
# Run source detection on the combined mosaic `{root}-ir_dr[cz]_sci.fits` and generates a catalog and segmentation image.  
# 
# Then perform simple matched-aperture photometry on the different available filter mosaics (in this case F098M and F140W from the direct imaging).  In principle the template fitting code shown below can incorporate this photometric information, though that's not currently done by default.

# In[30]:


## Run SEP (~SExtractor clone) catalog on the "ir" combined image
## and generate a photometric catalog with aperture photometry in all available bands
if not os.path.exists('{0}_phot.fits'.format(root)):
    get_background=False # SExtractor background subtraction
    tab = auto_script.multiband_catalog(field_root=root, threshold=1.8,
                                        detection_background=get_background,
                                        photometry_background=get_background) 
    
files = glob.glob('{0}-ir*'.format(root)) + glob.glob('*phot*fits')
for file in files:
    print(file)
    
phot = utils.GTable.gread('{0}_phot.fits'.format(root))
print('{0}Metadata{0}'.format('\n'+'='*20+'\n'))
for k in phot.meta:
    print('{0}:\t{1}'.format(k, phot.meta[k]))


# In[31]:


phot[:2].show_in_notebook()


# ## Building the grism exposure container: `multifit.GroupFLT` ##
# 
# With the preprocessing done, we can now start on the analysis of the spectra.  `Grizli` is built around low-level tools for modeling and analyzing each individual grism exposure individually.  Though once multiple exposures are available (e.g., exposures within a visit or separate visits with different grisms and/or orients) the information from each can be combined for analyzing the spectrum of a given object.  A benefit of the exposure-level processing is that all of the model-to-data comparisons (i.e. chi-squared) are done in the space of the original detector pixels, with their well-calibrated and well-understood noise properties.
# 
# The `GroupFLT` class provides a container for processing multiple FLT exposures simultanously.
# 
# ### Inputs ###
# * `grism_files` = list of grism exposure filenames
# * `direct_files` = (optional) list of direct exposure filenames
# * `ref_file` = (optional) reference direct image mosaic (one or the other of `ref_file` or `direct_files` should be specified.)
# * `seg_file`, `catalog` = segmentation image and catalog, usually generated with SExtractor
# * `cpu_count` = set to > 0 for parallel processing
# * `pad` parameter (default=200 pixels).  If set, then add padding around the FLTs to enable modeling of objects that would fall off of the direct image but that still disperse spectra onto the grism exposure (assuming they fall in the `ref_file` and `seg_file` mosaics).
# 
# The contents of the `grism_files` list can contain pretty much anything, with the result limited by memory / cpu power.  For example, you can provide a list of **all** 112 of the 3D-HST G141 exposures in the COSMOS field (4 exposures x 28 pointings), along with the field mosaic and segmentation images.  This example is actually fairly easy to process as individual objects will fall in typically 4, perhaps 8 individual exposures in some overlap areas.  Another example is a list of exposures from multiple instruments / grisms / orients of a single field, thought the subsequent fits can be slow if an object has spectra in many individual exposures.
# 
# Reference images are blotted to the distorted exposure frame with `AstroDrizzle.ablot`.   Messages go by, as below, when you load the `GroupFLT` object talking about "cutouts" because the script tries to make smaller cutouts of large reference images to speed up the blot processing.
# 
# **NB** Seems to often have memory leak problems if `seg_file` isn't significantly larger than the footprint of a given FLT file.  Drizzle `blot` segfaults out but the script just hangs since the multiprocessing threads don't get the message.
# 
# ### Flat continuum model ###
# 
# Once the `GroupFLT` object is initialized, compute a first-pass model of the full detector field of view assuming simple linear continua for all objects in the field (brighter than `mag_limit`).  By default this assumes a somewhat blue continuum suitable for the Rayleigh-Jeans tail of low-redshift objects.  It's far from perfect but the initial model does provide a good way of at least identifying which pixels of a given object could be contaminated by neighbors, even if the quantitative model is not precise.
# 
# ### Refined polynomial continuum model ###
# 
# After computing the simple continuum model, refine the model spectra for brighter objects using higher-order polynomials fit directly to the spectra themselves.  The `refine_list` method marches through objects starting with the brightest and fits a polynomial of order `poly_order` to the observed spectrum after subtracting off the model for contaminants.  Note that if the list of grism exposures contained multiple orientations covering a single field, this fit can be well constrained even in the presence of contamination.
# 
# The `grism_prep` script iterates on the refined polynomial model `refine_niter` times.
# 
# ### Save state ###
# 
# You can optionally dump saved data (i.e., `grp.save_full_data()`) for fast restart and avoid recomputing the contamination models, for example in a new Python session.  This can be done at any time after you've made changes to the GroupFLT data that you'd like to store for later.  The `grism_prep` script does this automatically.

# In[32]:


files = glob.glob('*GrismFLT.fits')
if len(files) == 0:
    ### Grism contamination model
    os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

    # Which filter to use as direct image?  Will try in order of the list until a match is found.
    gris_ref = {'G141': ['F140W', 'F160W'], 
                'G102': ['F105W', 'F098M', 'F110W']}

    x = auto_script.grism_prep(field_root=root, refine_niter=3,
                                 gris_ref_filters=gris_ref)

    grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), 
                            catalog='{0}-ir.cat.fits'.format(root), 
                            cpu_count=-1, sci_extn=1, pad=256)
    
else:
    grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), 
                            catalog='{0}-ir.cat.fits'.format(root), 
                            cpu_count=-1, sci_extn=1, pad=256)


# ### The final contamination model ###

# In[33]:


# Show the results of the contamination model
### Show FLT residuals



# ### Parameters for object fitting
# ### Read in LDP and match

# In[35]:


## USE THIS ONE
h = open('/data2/jrcooper/notebooks/reduction/EDisCS/gprior/j122816m1132s.txt', 'r')
lines = h.readlines()[1:]
h.close()   
z_LDP       = []  
id_HST      = []   
ra_HST      = []
dec_HST     = [] 

for line in lines: 
    a = line.split()     
    z_LDP.append(float(a[3]))
    ra_HST.append(float(a[1]))
    dec_HST.append(float(a[2]))



z_LDP     = np.array(z_LDP)
ra_HST    = np.array(ra_HST)
dec_HST   = np.array(dec_HST)
id_HST     = np.array(id_HST)   


# In[36]:


#### Store fit parameters to `fit_args.npy` for batch-mode processing

# Drizzle parameters for line maps
pline = auto_script.DITHERED_PLINE
#print(pline) 

sig = 0.007 
z = np.arange(0,2,.001)
for i in z_LDP:
    p_z = np.exp(-(z - i)**2/(2*sig**2))/((2*np.pi)**0.5/sig)

spec_prior = [z,p_z]
tuple
a = tuple(spec_prior)

# Generate the parameter dictionary
args = auto_script.generate_fit_params(field_root=root, prior=a, 
                                       MW_EBV=tabs[0].meta['MW_EBV'], 
                                pline=pline, fit_only_beams=True, run_fit=True, poly_order=7, 
                                fsps=True, sys_err = 0.03, fcontam=0.2, zr=[0.05, 3.4], 
                                save_file='fit_args.npy')


# ### Field PSF file ### 
# Make an average effective PSF for each available IR filter by evaluating the field-dependent PSF across the final mosaic and drizzling to a common output.  Also make an extension with a PSF on the pixel grid of the drizzled line map parameters generated above (`pline`).  Each PSF is generated with the native pixel grid and 2/4x oversampling for use with, e.g., [GALFIT](https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html).
# 
# *NB* There is currently no ePSF for F098M, so F105W is used (http://www.stsci.edu/~jayander/STDPSFs/WFC3IR/).
# 

# In[37]:


# Make PSF file
if not os.path.exists('{0}-f105w_psf.fits'.format(root)):
    auto_script.field_psf(root=root, HOME_PATH=HOME_PATH)


# In[38]:


# Show the PSFs
print(glob.glob('*psf.fits'),'\n')

im = pyfits.open('{0}-f105w_psf.fits'.format(root))
im.info()

fig = plt.figure(figsize=[6,6])
for i, ext in enumerate([1,2,4,5]):
    ax = fig.add_subplot(2,2,i+1)
    ax.imshow(np.log10(im[ext].data))
    
    h = im[ext].header
    label = '{0}-{1}  {2}"/pix \npixf={3} kernel={4} \nfilter={5}'.format(h['EXTNAME'],
                     h['EXTVER'], h['PSCALE'], h['PIXFRAC'], h['KERNEL'], h['FILTER'])
    
    ax.text(0.08, 0.95, label, ha='left', va='top', 
            transform=ax.transAxes, size=8, color='k',
            bbox={'fc':'w'})
    
    ax.set_xticklabels([]); ax.set_yticklabels([])
    
    sh = im[ext].data.shape
    ax.set_xlim(sh[1]/2-40, sh[1]/2+40)
    ax.set_ylim(sh[0]/2-40, sh[0]/2+40)
    
fig.tight_layout(pad=1)


# ## Extract and fit individual spectra ##

# In[39]:


os.chdir('../Extractions')


# In[40]:


### Find IDs of specific objects to extract
import grizli.pipeline
from grizli.pipeline import auto_script
import astropy.units as u
tab = utils.GTable()
#tab['ra'] = [174.4118633,174.3686397,174.3797582,174.3845889,174.397861,174.412518,174.3670438,174.3742016,174.3804661,174.4072525,174.4070533,174.3933064,174.3931596,174.3691858,174.3688672,174.3690021,174.3688063,174.3672425,174.3665378,174.3921713,174.3933575,174.4087276,174.4123249,174.3767788,174.3748348,174.3787987,174.3708528,174.3812335,174.3814446,174.382614,174.3791532,174.385079,174.3849722,174.3773103,174.3779308,174.3717816,174.3848934,174.3771574,174.3727826,174.3734838,174.3913213,174.3869972,174.3864554,174.3866813,174.3799761,174.3879,174.3961038,174.3945939,174.3946899,174.3980747,174.3941714,174.3810396,174.3811199,174.3759067,174.3738911,174.3804884,174.3793689,174.3769608,174.3739911,174.3941685,174.4011006,174.3778083,174.3928164,174.3804512,174.3851616,174.380016,174.3801094,174.3965471,174.3884755,174.3910637,174.387964,174.3913192,174.3914445,174.3802608,174.3800891,174.3800327,174.390005,174.3766207,174.3861459,174.3751421,174.4029965,174.3937154,174.3938229,174.4006638,174.400634,174.4004205,174.3959245,174.3834155,174.410023,174.4110051,174.3941949,174.4079271,174.3832748,174.3723437,174.3752115,174.3881324,174.40738,174.4075371,174.4038001,174.3776531,174.3985774,174.4051396,174.4012677,174.3806377,174.3716552,174.3846171,174.3850321,174.3871291,174.3974057,174.3707338,174.4001625,174.4103907,174.4104152,174.370494,174.3832325,174.3716295,174.3823567,174.3849623,174.4010585,174.3748725,174.4069045,174.3850529,174.4107016,174.3988903,174.3771158,174.4048342,174.373967,174.4067188,174.4034601,174.4040355,174.4038954,174.3922065,174.3925352,174.4026382,174.3995775,174.3831676,174.3978151,174.394407,174.3943988,174.3804665,174.3799555,174.4073365,174.3810505,174.3989021,174.4118035,174.3862916,174.3950565,174.3815119,174.4114332,174.4089518,174.4086337,174.3783956,174.3835116,174.3833891,174.3772125,174.3767481,174.3763366,174.3686133,174.3846834,174.3765115,174.3710509,174.3811286,174.3810785,174.3878796,174.3714472,174.3696215,174.3837304,174.3799637,174.3975855,174.408378,174.389328,174.4031666,174.3689918,174.3821339,174.3731581,174.3684657,174.3997755,174.3966998,174.4092418,174.3988337,174.4026826,174.399456,174.3943523,174.4099742,174.4053917,174.381934,174.3701956,174.3878967,174.403905,174.411098,174.3935612,174.4090161,174.387054,174.3865922,174.3829451,174.3755649,174.4105957,174.3842971,174.3809018,174.409077,174.408885,174.3994039,174.3691348,174.3796697,174.4095749,174.4049898,174.378185,174.387129,174.4055909,174.4055366,174.4057687,174.4054303,174.4003571,174.3982444,174.4062903,174.4061453,174.3997543,174.3820729,174.3840378,174.3682332,174.3686231,174.3949571,174.3994276,174.3875546,174.3964459,174.4050611,174.394733,174.3806533,174.3815379,174.3955579,174.3894354,174.401443,174.3935624,174.3934467,174.4084085,174.401406,174.3988417,174.3707744,174.4083716,174.4092932,174.4034055,174.3907952,174.3682603,174.3968084,174.395399,174.3981921,174.3752158,174.3711835,174.3986423,174.3842199,174.4028942,174.3918513,174.3917949,174.3996218,174.4063234,174.3943757,174.3907661,174.3867227,174.4020477,174.4020128,174.3996745,174.3861347,174.3837594,174.4056763,174.3677452,174.3945288,174.4021971,174.4080006,174.3782513,174.3986636,174.3678685,174.3740759,174.3721,174.3722116,174.3823471,174.3954012,174.3750514,174.3992077,174.4014894,174.3749131,174.407115,174.3823573,174.3673895,174.3785874,174.3838518,174.400227,174.4071072,174.4071757,174.3945508,174.3776495,174.3726168,174.3786074,174.4013678,174.3726654,174.4058475,174.3865453,174.3724778,174.3747123,174.3748985,174.3801947,174.4056129,174.399787,174.3948029,174.3943262,174.3814047,174.4094231,174.3891988,174.3785227,174.4020965,174.366548,174.4030109,174.4029031,174.3710936,174.3701846,174.399441,174.3849507,174.3997734,174.3855422,174.385701,174.3940045,174.3774031,174.404311,174.3888625,174.368653,174.392602,174.386554,174.3951195,174.3862442,174.3884426,174.3933544,174.4037673,174.3962347,174.3957018,174.3962323,174.3742381,174.3999549,174.3997852,174.3847455,174.3748441,174.3749965,174.3827494,174.36688,174.3972047,174.4012877,174.4033545,174.4032237,174.3922029,174.3793299,174.3700605,174.3941671,174.3776026,174.3996319,174.3814074,174.4024752,174.4025639,174.3715811,174.3735299,174.380875,174.3824583,174.377477,174.3731738,174.373008,174.3702185,174.3953243,174.3863515,174.3702485,174.4002528,174.4029622,174.3655957,174.3892633,174.39352,174.3948536,174.393928,174.377719,174.3672721,174.406376,174.3657231,174.3972316,174.3684152,174.3808368,174.3739879,174.397792,174.399329,174.4060634,174.3955529,174.389342,174.3835584,174.4032463,174.3726363,174.3902532,174.3832384,174.3968654,174.4053413,174.3734288,174.3908544,174.3691144,174.3745796,174.3744327,174.3707195,174.3909923,174.3763947,174.3730874,174.3707268,174.4030591,174.3910123,174.3949013,174.4014471,174.3801271,174.4053426,174.3658806,174.4046002,174.4019161,174.4015301,174.3906148,174.3731218,174.3860184,174.3821069,174.3916208,174.3915034,174.3736736,174.4044422,174.3705266,174.3983853,174.397938,174.3769437,174.3816981,174.3854264,174.3685865,174.3764316,174.3735563,174.4045562,174.4047738,174.3661968,174.3821705,174.4051944,174.383714,174.3684348,174.404838,174.3959258,174.3956777,174.382061,174.3995164,174.3766376,174.4025186,174.3805095,174.3953532,174.3686929,174.3688924,174.3685074,174.3835799,174.4012345,174.3854665,174.366571,174.3984712,174.3736616,174.3907737,174.384933,174.3682905,174.3673654,174.3947067,174.3956914,174.4015167,174.366419,174.3663184,174.3677471,174.3913106,174.3816995,174.3816084,174.3918781,174.3677641,174.4021059,174.3903793,174.3807377,174.3831399,174.4066423,174.3762605,174.3871168,174.3928499,174.3785797,174.4050567,174.4037366,174.3746231,174.395452,174.3842919,174.4009807,174.4007248,174.3803011,174.3877482,174.3841823,174.392275,174.3828127,174.404194,174.3887049,174.3810874,174.3957987,174.3999629,174.3999316,174.3852252,174.3947081,174.393336,174.3962873,174.4037769,174.4025584,174.3950435,174.3972366,174.3896022,174.4014372,174.3985801,174.4042522,174.4056772,174.4008711,174.4036353,174.4023552]
#tab['dec'] = [-11.42853799,-11.43392584,-11.44203982,-11.44097895,-11.43806378,-11.43475228,-11.40947934,-11.40777854,-11.40653533,-11.41045374,-11.40955852,-11.43885947,-11.43889343,-11.40906908,-11.40907182,-11.4347663,-11.43440486,-11.42673946,-11.42402156,-11.43908326,-11.4387369,-11.43550106,-11.43466296,-11.4415287,-11.44184637,-11.44110515,-11.44092653,-11.44075424,-11.44073077,-11.44090179,-11.44046386,-11.44052153,-11.44043305,-11.44010887,-11.43995076,-11.43966209,-11.43957531,-11.4391517,-11.43900826,-11.43880573,-11.43874973,-11.43845038,-11.43867664,-11.43854369,-11.43819859,-11.43797968,-11.4376107,-11.43791397,-11.43765249,-11.43769736,-11.43748824,-11.43737302,-11.4373696,-11.43733843,-11.43729457,-11.43737927,-11.43731216,-11.43716304,-11.43700523,-11.43715061,-11.43695245,-11.43738873,-11.43667961,-11.43660114,-11.43664056,-11.4368921,-11.43669941,-11.43640085,-11.43636224,-11.43614048,-11.43615751,-11.43621187,-11.4360709,-11.43625207,-11.43621388,-11.43604582,-11.43592263,-11.43581439,-11.43547686,-11.4354173,-11.43523482,-11.43517732,-11.43513747,-11.43511685,-11.43503007,-11.43502985,-11.43511864,-11.43511806,-11.43490389,-11.43487456,-11.43484514,-11.43562152,-11.43480915,-11.43491464,-11.43477839,-11.43473346,-11.435125,-11.43478885,-11.43450863,-11.4344837,-11.43407238,-11.43416257,-11.43422246,-11.43398363,-11.4339249,-11.43405881,-11.43413046,-11.43387524,-11.43378928,-11.43378818,-11.43354809,-11.43384645,-11.43361901,-11.43329336,-11.43371806,-11.43326885,-11.43317077,-11.43311573,-11.43309383,-11.43334678,-11.43294878,-11.43297176,-11.43254428,-11.43249794,-11.4324047,-11.43238153,-11.43246462,-11.43243397,-11.4324398,-11.43246459,-11.43238497,-11.43315044,-11.43261653,-11.43210802,-11.43205931,-11.43208309,-11.43196333,-11.43262748,-11.43217416,-11.43225441,-11.43230031,-11.4315085,-11.43148172,-11.43134974,-11.43134326,-11.43122554,-11.43114979,-11.43111243,-11.43107744,-11.43195058,-11.43130636,-11.43085946,-11.4309923,-11.43080579,-11.43107831,-11.43115216,-11.43162236,-11.43083445,-11.43064685,-11.43052335,-11.43049765,-11.43052102,-11.43044758,-11.4304767,-11.43026676,-11.43032658,-11.43025434,-11.43036853,-11.43009585,-11.4299374,-11.42997225,-11.42979159,-11.42926305,-11.42930552,-11.42928875,-11.42927311,-11.42913779,-11.42911706,-11.42920944,-11.42908388,-11.42906351,-11.42903199,-11.42903966,-11.428905,-11.42885695,-11.42880471,-11.42879791,-11.42879156,-11.4287385,-11.4290915,-11.42882846,-11.4286187,-11.42860604,-11.42850556,-11.42853758,-11.42846891,-11.42839968,-11.42833504,-11.42823358,-11.42813814,-11.4279825,-11.42794405,-11.42793854,-11.42789301,-11.42790229,-11.42788122,-11.42786879,-11.42785784,-11.42781333,-11.42760363,-11.42768704,-11.4279682,-11.42738644,-11.42729563,-11.42734539,-11.42734684,-11.42725568,-11.4272588,-11.42723385,-11.42708004,-11.42701326,-11.42675426,-11.4272413,-11.42673909,-11.42669072,-11.42663224,-11.42653873,-11.42645318,-11.42641739,-11.42655658,-11.42630766,-11.42620436,-11.42640751,-11.42625184,-11.42613881,-11.42602471,-11.42594957,-11.42595832,-11.42584453,-11.42578177,-11.42569227,-11.42563155,-11.42562447,-11.42506475,-11.42467519,-11.42462967,-11.42458801,-11.42452296,-11.424535,-11.42417356,-11.42398218,-11.42404442,-11.42395455,-11.4239113,-11.42380833,-11.4241093,-11.42374028,-11.42372406,-11.4236061,-11.42353652,-11.42352977,-11.42356248,-11.42347299,-11.42346748,-11.42341684,-11.42333948,-11.42322254,-11.42321213,-11.42338193,-11.4231503,-11.42315317,-11.42316601,-11.42307617,-11.42304559,-11.42307567,-11.42301023,-11.42358478,-11.42286573,-11.42277898,-11.42269527,-11.42277888,-11.42261934,-11.42257668,-11.4224877,-11.42248283,-11.42266092,-11.42249751,-11.42238369,-11.4223102,-11.42228857,-11.42223818,-11.4222577,-11.42221041,-11.42197042,-11.42192727,-11.42192789,-11.42188319,-11.4221875,-11.42196477,-11.42174339,-11.42169342,-11.42159139,-11.42159263,-11.42154861,-11.42155793,-11.42149116,-11.42133729,-11.42123386,-11.42122765,-11.42128559,-11.42134311,-11.42129676,-11.42118796,-11.42106807,-11.42105894,-11.42097638,-11.42090431,-11.42085083,-11.42083516,-11.42108404,-11.42086033,-11.42080707,-11.42080744,-11.42071458,-11.42064634,-11.42068774,-11.42052511,-11.42045126,-11.42028535,-11.42025722,-11.42032373,-11.42001829,-11.42040404,-11.42030859,-11.41992769,-11.41964469,-11.4196182,-11.41981048,-11.41922171,-11.41919473,-11.41918333,-11.41914083,-11.41907501,-11.41893541,-11.4188827,-11.41890119,-11.4188498,-11.41853853,-11.41854232,-11.4181411,-11.41815863,-11.41800258,-11.41799289,-11.41834506,-11.4180769,-11.41792495,-11.41782528,-11.41823586,-11.41768962,-11.41773111,-11.41771763,-11.41746561,-11.41740475,-11.41737048,-11.4174037,-11.4172372,-11.41729615,-11.41717035,-11.41710237,-11.41709022,-11.41702467,-11.41689284,-11.41685202,-11.4167902,-11.41675159,-11.41694998,-11.41661692,-11.41660054,-11.41650823,-11.41651522,-11.41624324,-11.41612911,-11.41589484,-11.41570574,-11.41570162,-11.4157157,-11.41566355,-11.41567259,-11.41559252,-11.41568055,-11.41575505,-11.41545817,-11.41544078,-11.41539232,-11.41533708,-11.41510099,-11.41503168,-11.41501003,-11.41507477,-11.41482434,-11.41484798,-11.41462203,-11.41462011,-11.41464731,-11.41454381,-11.41450683,-11.41449022,-11.41448697,-11.41442341,-11.41438362,-11.41437638,-11.41431681,-11.41419201,-11.41402041,-11.41402811,-11.41405514,-11.4138771,-11.41386694,-11.41375663,-11.41375997,-11.41364712,-11.41350653,-11.413318,-11.41314819,-11.41305758,-11.41288839,-11.4127817,-11.41266544,-11.41259104,-11.41246816,-11.41291447,-11.41268732,-11.41237211,-11.41235695,-11.41226611,-11.41216348,-11.4121622,-11.41180646,-11.41186404,-11.41175891,-11.41189794,-11.41163934,-11.41149071,-11.41140638,-11.41135733,-11.41131546,-11.4116813,-11.41154078,-11.41161017,-11.41123261,-11.41090604,-11.41075214,-11.41071288,-11.41066829,-11.41067879,-11.41091285,-11.41052618,-11.41051927,-11.41046357,-11.41045116,-11.41039976,-11.41032549,-11.41039541,-11.41028855,-11.41025512,-11.41020973,-11.41010301,-11.41000598,-11.40994451,-11.4098655,-11.40970717,-11.40974273,-11.40965623,-11.40955383,-11.40941978,-11.40938211,-11.40929036,-11.40920089,-11.40934505,-11.40864362,-11.40844014,-11.40846545,-11.40839585,-11.40831762,-11.40849842,-11.40816308,-11.40791683,-11.40792574,-11.40781022,-11.40795797,-11.40788815,-11.40735788,-11.40744121,-11.40703858,-11.40682265,-11.40669668,-11.40679066,-11.40643416,-11.40644136,-11.40572145,-11.40547421,-11.40520698,-11.40522151,-11.40546231,-11.40514909,-11.40515372,-11.40511544,-11.40476582,-11.40521459,-11.4055198,-11.40445166,-11.40359569,-11.40297454]
tab['ra'] = ra_HST
tab['dec'] = dec_HST
idx, dr = grp.catalog.match_to_catalog_sky(tab)
source_ids = grp.catalog['NUMBER'][idx]
tab['id'] = source_ids
tab['dr'] = dr.to(u.mas)
tab['dr'].format='.1f'
tab.show_in_notebook()


# ### Extract 2D spectra "beams" ###
# 
# The `GroupFLT` object contains the entire exposure information, from which we can make cutouts of spectra for individual objects with the `get_beams` method.  These cutouts are more managable and portable than the entire exposures, though currently the processing does work in the paradigm of having a static contamination model for a given object.  
# 
# In pipeline mode, the function below is called with `ids=[], maglim=[mag_min, mag_max]` and all objects in the reference catalog with `mag_min < MAG_AUTO < mag_max` are extracted.  The redshift fits are performed if `run_fit=True`.

# In[41]:


for id_i in source_ids:
    auto_script.extract(field_root=root, prior=a, MW_EBV=tabs[0].meta['MW_EBV'], 
                                pline=pline, fit_only_beams=True, run_fit=True, poly_order=7, 
                           grp=grp, diff=True)


# In[49]:

print('all done!')
#id=source_ids[0]
#for id_i in source_ids:
#    auto_script.extract(field_root=root, ids=[id_i], MW_EBV=tabs[0].meta['MW_EBV'], 
  #                  pline=pline, run_fit=True, grp=grp, diff=True)


# ### 2D spectra ###
# 
# The spectral extraction produces two versions of the extracted 2D spectra:
# 
# * `{root}_{id:05d}.beams.fits` : Multi-extension FITS file with sets of extensions for 2D cutouts **from each individual grism exposure**.  Fitting in this space is most robust as the grism dispersion is defined in the "FLT" coordinates and the model comparison is done directly on un-resampled image pixels with relatively well-understood noise properties.
#     
#     
# * `{root}_{id:05d}.stack.fits` : Multi-extension FITS file with extension with combinations of all exposures in a given grism & position angle.  The fitting tools can be used with these products as well, where the fits are much faster as 2D models at each trial redshift are produced for `N_PA x N_grism` combinations, where often `N_PA x N_grism << N_exposure`.   The fits are evaluated in the resampled drizzled pixel space, and they are often less robust than fits to the full "beams" spectra, particularly at low S/N.
#     
#     The `{root}_{id:05d}.stack.png` files, shown below, are often useful for visual inspection of the 2D products.  Note that the bottom panel of the `stack.png` files is the drizzled combination of *all* PAs for a given grism, and with a polynomial continuum model subtracted if `diff=True` in the extraction script above.
#     

# In[40]:


#Image(filename='j113812m1134_00362.stack.png')
#Image(filename='{0}_{1:05d}.stack.png'.format(root, id)) 


# In[41]:


# 1D spectrum with polynomial model
#Image(filename='j113812m1134_00362.1D.png')
#Image(filename='{0}_{1:05d}.1D.png'.format(root, id)) 


# ### Redshift fit ###
# 
# The redshift fit is performed in the following steps:
# 
# * On a coarse redshift grid (dz/1+z ~ 0.005) fit continuum templates along with **line complex** templates for a) [OII]+[NeIII], b) [OIII]+Hbeta, and c) Halpha+[SII]+weaker red lines.  These line complexes have fixed line ratios but are useful for breaking redshift degeneracies as these lines do, usually, come in groups.  Leaving all line strengths free would allow for perfect degeneracy between, e.g., Halpha and [OII] (assuming no significant continuum features).
# 
# * Find peaks (minima) in the chi-squared on the coarse grid and zoom in on them now allowing for more freedom in the indifidual line strengths, as well as fitting on a fine redshift grid sufficient to resolve the best redshift.
# 
# **NB** Continuum templates are needed in the directory `${GRIZLI}/templates`.  The template names are currently hard-coded in [multifit.py](https://github.com/gbrammer/grizli/blob/master/grizli/multifit.py) and the easiest way to make them available is to symlink them from the `data/templates` directory that accompanies the `grizli` code distribution:
# 
# 
# ### Emission line maps ###
# 
# Once we've computed the full continuum + line model, we can create 2D *drizzled* maps at any desired output wavelength, for example to make emission line maps.  This makes use of the WCS information in the individual grism FLT exposures and the outputs can have any desired WCS (e.g., pixel scale & dimensions) and can be used to compare directly to imaging data.
# 
# The emission line maps are generated by subtracting the best-fit continuum model, assuming that the direct image is representative of the continuum morphology.  This should be a reasonable assumption for objects other than, perhaps, those with extreme line equivalent widths.
# 

# In[44]:


# Fit it.  The "run_all_parallel" function defaults to all of the parameters set in 'fit_args.npy'
#for id_i in source_ids:
#    fitting.run_all_parallel(id_i)


# ### Fit products ###
# 
# A number of files are produced that contain the results of the redshift fit.  The [`NewSpectrumFits.ipynb`](https://github.com/gbrammer/grizli/blob/master/examples/NewSpectrumFits.ipynb) notebook describes how to interact with these products in some greater detail.  

# In[154]:


#files = glob.glob('*{0:05d}*'.format(id))
#for file in files:
#    print(file)


# In[155]:


#for file in files:
#    if not file.endswith('.fits'):
 #       continue
        
 #   im = pyfits.open(file)
 #   print('\n\n{1}\n{0}\n{1}\n\n'.format(file, '='*len(file)))
 #   im.info()


# ### Continuum-dominated spectra ###
# 
# The object below is the dominated by strong Balmer break and absorption lines (see [van Dokkum & Brammer 2010](http://adsabs.harvard.edu/abs/2010ApJ...718L..73V)).  The redshift fit and spectral constraints are precise even without any supporting photometric data.

# In[156]:


# Continuum source
#id=source_ids[0:512]
#id=np.array(id)
#auto_script.extract(field_root=root, ids=id, MW_EBV=tabs[0].meta['MW_EBV'], 
          #          pline=pline, run_fit=False, grp=grp, diff=True)


# In[142]:


# Stacked 2D spectrum
#Image(filename='{0}_{1:05d}.stack.png'.format(root, id)) 


# In[143]:


# 1D spectrum with polynomial model
#Image(filename='{0}_{1:05d}.1D.png'.format(root, id)) 


# In[45]:


## Run the fit
#for id_i in source_ids:
  #  fitting.run_all_parallel(id_i)


# ### Fit grism with photometry
# 
# Another option is fitting the grism spectra along with ancillary photometry, described here: [Fit-with-Photometry.ipynb](https://github.com/gbrammer/grizli/blob/master/examples/Fit-with-Photometry.ipynb).
