from enthought.mayavi.core.registry import registry
from enthought.mayavi.core.metadata import SourceMetadata
from enthought.mayavi.core.pipeline_info import PipelineInfo


def treader(fname, engine):
    from enthought.tvtk.api import tvtk
    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
    from enthought.mayavi.sources.api import ArraySource
    # Do your own reader stuff here, I'm just reading a VTK file with a
    # different extension here.
    
    src = ArraySource()
    src.scalar_name=fname
    import os.path
    src.name=os.path.split(fname)[-1]
    return src

treader_info = SourceMetadata(
    id            = "Sp4Array File Reader",
    factory = 't_reader.treader',
    tooltip       = "Load a Sp4 file",
    desc   = "Load a Sp4 file",
    help   = "Load a Sp4 file",
    menu_name        = "&Sp4 file",
    extensions = ['sp4','SP4'],
    wildcard = 'Sp4 files (*.Sp4)|*.Sp4',
    output_info = PipelineInfo(datasets=['unstructured_grid'],
                               attribute_types=['any'],
                               attributes=['any'])
)
# Inject this information in the mayavi registry
registry.sources.append(treader_info)


if __name__=='__main__':
	source=tSource()
	source.configure_traits()
