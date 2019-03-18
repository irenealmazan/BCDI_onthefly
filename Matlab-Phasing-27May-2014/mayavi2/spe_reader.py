from enthought.mayavi.core.registry import registry
from enthought.mayavi.core.metadata import SourceMetadata
from enthought.mayavi.core.pipeline_info import PipelineInfo


def spe_reader(fname, engine):
    """Reader for .zzz files.
    
    Parameters:
    -----------
    fname -- Filename to be read.
    engine -- The engine the source will be associated with.
    """
    from enthought.tvtk.api import tvtk
    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
    # Do your own reader stuff here, I'm just reading a VTK file with a
    # different extension here.
    import SpeFile as spe
    r = spe.SpeFile(fname)

    from enthought.mayavi.sources.api import ArraySource
    src = ArraySource()
    src.scalar_data=r.GetNumArray()
    src.scalar_name=fname
    import os.path
    src.name=os.path.split(fname)[-1]
    return src

spe_reader_info = SourceMetadata(
    id            = "Winview File Reader",
    factory = 'spe_reader.spe_reader',
    tooltip       = "Load a SPE file",
    desc   = "Load a SPE file",
    help   = "Load a SPE file",
    menu_name        = "&SPE file",
    extensions = ['spe','SPE'],
    wildcard = 'SPE files (*.SPE)|*.SPE',
    output_info = PipelineInfo(datasets=['image_data'],
                               attribute_types=['any'],
                               attributes=['any'])
)
# Inject this information in the mayavi registry
registry.sources.append(spe_reader_info)
