from sage.misc.decorators import options
@options(plot_points=20)
def plot_vector_field_on_curve( (xf, yf), (x, y), range, **options ):
    r"""Plot values of a vector-values function along points of a curve
    in the plane.

    Note this function doesn't plot the curve itself."""
    from sage.plot.all import Graphics
    from sage.misc.misc import xsrange
    from sage.plot.plot_field import PlotField
    from sage.plot.misc import setup_for_eval_on_grid
    zz, rangez = setup_for_eval_on_grid( (x, y, xf, yf), [ range ], options['plot_points'] )
    #print 'setup: ', zz, rangez
    x, y, xf, yf = zz
    xpos_array, ypos_array, xvec_array, yvec_array = [],[],[],[]
    for t in xsrange( *rangez[0], include_endpoint=True ):
       xpos_array.append( x(t) )
       ypos_array.append( y(t) )
       xvec_array.append( xf(t) )
       yvec_array.append( yf(t) )
    import numpy
    xvec_array = numpy.ma.masked_invalid(numpy.array(xvec_array, dtype=float))
    yvec_array = numpy.ma.masked_invalid(numpy.array(yvec_array, dtype=float))
    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    g.add_primitive(PlotField(xpos_array, ypos_array, xvec_array, yvec_array, options))
    return g

