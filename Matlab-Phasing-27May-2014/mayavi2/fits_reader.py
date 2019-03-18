# Based on the Mayavi2 zzz_reader example, by James Turner and
# Prabhu Ramachandran, SciPy 2008.

"""This is a simple reader factory to register PyFITS with mayavi as
   a file reader.

To use this:

    - put this in ~/.mayavi2/
    - then import this module in your ~/.mayavi2/user_mayavi.py.

thats it.

What you should get:

    - Options to open .fits files from the file->open menu.
    - Open .fits files via right click.
    - Open .fits files from the engine or mlab (via open)
    - do mayavi2 -d foo.fits.

"""

from enthought.mayavi.core.registry import registry
from enthought.mayavi.core.metadata import SourceMetadata
from enthought.mayavi.core.pipeline_info import PipelineInfo

def FITS_reader(fname, engine):
    """Reader for .fits files.
    
    Parameters:
    -----------

    fname -- Filename to be read.

    engine -- The engine the source will be associated with.
    """
    from enthought.tvtk.api import tvtk
    from enthought.mayavi.sources.array_source import ArraySource

    # Import some stuff needed to deal with dialogue boxes:
#    from enthought.pyface.api import FileDialog, OK
#    from enthought.pyface.action.api import Action
#    from enthought.mayavi.core.common import error

    # Import PyFITS if available, otherwise generate an error:
#    try:
#        import pyfits
#    except:
#        print "ERROR: couldn't import PyFITS module"
#        print "  (http://www.stsci.edu/resources/software_hardware/pyfits)"
#        # error("", parent) ?? How do we open the dialogue box?
#        return
    import pyfits


          
    # Open the FITS file with PyFITS:
    hdulist = pyfits.open(fname)  # What do we do if it doesn't exist?

    # How many Image extensions does the FITS file have? If it's simple
    # FITS or has a single image HDU, load that, otherwise we have to pop
    # up a little GUI window to select which extension to load:

    # Look in ETS_3.0.0/Mayavi_3.0.0/enthought/mayavi/action/sources.py

    # Currently just open extension 1
    data = hdulist[1].data.astype('d')
    src = ArraySource(scalar_data=data)
    return src   # If we return None, it won't do anything, but I now have
                 # a reference to the engine, so I can pop up a GUI to
                 # select a FITS extension and then later call ArraySource
                 # (or just do it above)
    # Can also subclass ArraySource (eg. metadata inside mayavi source,
    # small changes to FITS_reader_info below)

FITS_reader_info = SourceMetadata(
    id            = "FITSReader",
    factory = 'fits_reader.FITS_reader',
    tooltip       = "Load a FITS file",
    desc   = "Load a FITS file",
    help   = "Load a FITS file",
    menu_name        = "&FITS file",
    extensions = ['fits'],
    wildcard = 'FITS files (*.fits)|*.fits',
    output_info = PipelineInfo(datasets=['unstructured_grid'],
                               attribute_types=['any'],
                               attributes=['any'])
)
# Inject this information in the mayavi registry
registry.sources.append(FITS_reader_info)

if __name__ == '__main__':
    import sys
    print "*"*80
    print "ERROR: This script isn't supposed to be executed."
    print __doc__
    print "*"*80
    sys.exit(1)

