def user_outline():
    """A Factory function that creates a new module to add to the
    pipeline.  Note that the method safely does any mayavi imports
    inside avoiding any circular imports.
    """
    print "User Outline"
    from enthought.mayavi.modules.outline import Outline
    o = Outline(outline_mode='cornered', name='UserOutline')
    return o
