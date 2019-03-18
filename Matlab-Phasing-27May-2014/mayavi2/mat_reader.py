from enthought.mayavi.core.registry import registry
from enthought.mayavi.core.metadata import SourceMetadata
from enthought.mayavi.core.pipeline_info import PipelineInfo

def mat_reader(fname, engine):
    """Reader for .zzz files.
    
    Parameters:
    -----------
    fname -- Filename to be read.
    engine -- The engine the source will be associated with.
    """

    import os.path

    dir=os.path.dirname(fname)
    phparams=os.path.join(dir, "phasingparams.py")

    from JesseMatLabCoordSource import JesseMatLabCoordSource
    src = JesseMatLabCoordSource(file_name=fname, engine=engine)
    src.scalar_name=fname
    src.name=os.path.split(fname)[-1]
    return src

mat_reader_info = SourceMetadata(
    id            = "Matlab File Reader",
    factory = 'mat_reader.mat_reader',
    tooltip       = "Load a mat file",
    desc   = "Load a mat file",
    help   = "Load a mat file",
    menu_name        = "&mat file",
    extensions = ['rec','REC','mat'],
    wildcard = 'Rec files (*.rec)|*.Rec',
    output_info = PipelineInfo(datasets=['image_data','structured_grid'],
                               attribute_types=['any'],
                               attributes=['any'])
)
# Inject this information in the mayavi registry
registry.sources.append(mat_reader_info)


if __name__=='__main__':
	source=Sp4ArrayFileSource(file_name='/Users/rharder/PhasingProjects/diffmap-moyuAu/Au708-81-Imask400rs/BD-data.sp4')
	source.update()
	print source._scalar_data
	source.configure_traits()
