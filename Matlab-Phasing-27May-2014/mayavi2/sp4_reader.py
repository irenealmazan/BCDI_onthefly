from enthought.mayavi.core.registry import registry
from enthought.mayavi.core.metadata import SourceMetadata
from enthought.mayavi.core.pipeline_info import PipelineInfo

def sp4_reader(fname, engine):
    """Reader for .zzz files.
    
    Parameters:
    -----------
    fname -- Filename to be read.
    engine -- The engine the source will be associated with.
    """

    import os.path

    dir=os.path.dirname(fname)
    phparams=os.path.join(dir, "phasingparams.py")

    if os.path.isfile(phparams):
        from Sp4ArrayCoordSource import Sp4ArrayCoordSource
    	src = Sp4ArrayCoordSource(file_name=fname, engine=engine)
    	src.scalar_name=fname
    	src.name=os.path.split(fname)[-1]
    else:
	from Sp4ArraySource import Sp4ArraySource
	src = Sp4ArraySource(file_name=fname, engine=engine)
    return src

sp4_reader_info = SourceMetadata(
    id            = "Sp4Array File Reader",
    factory = 'sp4_reader.sp4_reader',
    tooltip       = "Load a Sp4 file",
    desc   = "Load a Sp4 file",
    help   = "Load a Sp4 file",
    menu_name        = "&Sp4 file",
    extensions = ['sp4','SP4'],
    wildcard = 'Sp4 files (*.Sp4)|*.Sp4',
    output_info = PipelineInfo(datasets=['image_data','structured_grid'],
                               attribute_types=['any'],
                               attributes=['any'])
)
# Inject this information in the mayavi registry
registry.sources.append(sp4_reader_info)


if __name__=='__main__':
	source=Sp4ArrayFileSource(file_name='/Users/rharder/PhasingProjects/diffmap-moyuAu/Au708-81-Imask400rs/BD-data.sp4')
	source.update()
	print source._scalar_data
	source.configure_traits()
