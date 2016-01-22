"""Dynamical systems classes"""
from sage.all import *
from string import join
from sage.misc.latex import _latex_file_
latex.add_to_preamble('\\usepackage{amsmath}')

# =============================================================================
# Core functions and classes for dynamical systems processing in Sage
# =============================================================================

# First helper functions and classes

class DynamicsException(Exception):
    def __init__( self, message, latex_str=None ):
        self._latex_str = latex_str
        Exception.__init__( self, message )

class EquilibriumDetectedException(DynamicsException):
    def __init__(self):
        Exception.__init__( self, 'Equilibrium detected' )

class UnboundedDynamicsException(DynamicsException):
    def __init__(self, message=None):
        if message is None:
            message = 'Dynamics has become unbounded'
        Exception.__init__( self, message )

class TrajectoryInterruptedException(DynamicsException):
    # this one is kind of meta. It says the integration of the dynamics
    # has been interrupted, and provides the partial trajectory together
    # with the reason why - and the reason why is another exception...
    def __init__( self, trajectory, reason ):
        self._trajectory = trajectory
        self._reason = reason
        Exception.__init__( self, 'Integration of dynamics interrupted' )

def xform_symbol(v, xform_str, xform_latex):
    """little utility function to do things like add a hat or similar thing
    to a Sage variable, e.g. make X_i into \hat{X}_i."""
    try:
        # is it a symbol?
        if not v.is_symbol():
            # no, it's a more complex expression
            raise ValueError( str(v) + ' is not a symbol' )
    except AttributeError:
        # it's a string
        v = SR.symbol(v)
    import re
    name = str(v)
    mn = re.match( '^([^^_]*)(.*?)$', name )
    if mn: name = xform_str( mn.group(1), mn.group(2) )
    #name = re.sub( '^([^_]*)(.*?)$', r'\g<1>'+flair+'\g<2>', str(v) )
    latex_name = latex(v)
    #latex_name = re.sub( r'^([^_]*)(.*?)$', r'\\'+flair+'{\g<1>}\g<2>', latex(v) )
    ml = re.match( r'^([^^_]*)(.*?)$', latex_name )
    latex_name = xform_latex( ml.group(1), ml.group(2) )
    return SR.symbol( name, latex_name=latex_name )

def add_flair(v, flair):
    def xform_str( base, subs ): return base + flair + subs
    def xform_latex( base, subs ): return r'\\' + flair + '{' + base + '}' + subs
    return xform_symbol(v, xform_str, xform_latex)

def hat(v):
    return add_flair(v, 'hat')

def dot(v):
    return add_flair(v, 'dot')

class Trajectory(SageObject):
    """A trajectory of a system of ordinary differential equations.

    This class remembers the ODE that made it, and the symbolic variables
    involved in the ODE, so that the trajectories of complex symbolic
    expressions can be evaluated and plotted."""
    def __init__(self, system, timeseries):
        """system: the ODESystem object that provided the dynamics
        timeseries: a list of dictionaries, each providing a set of
            variable->number mappings, including the time variable.
            or a list of lists/tuples of values, ordered as the state
            variables are, with time prepended to each.
        """
        self._system = system
        try:
            timeseries[0]( system._time_variable ) # check if it's a bindings
        except TypeError:
            # make that timeseries of lists into timeseries of bindings
            #print 'make timeseries from', timeseries
            timeseries = [
                self.point_to_bindings( point ) for point in timeseries
            ]
            #print 'get', timeseries
	except IndexError: # timeseries is empty
	    pass
        self._timeseries = timeseries
    def __repr__(self):
        return 'Trajectory(' + repr(self._timeseries) + ')'
    def point_to_bindings( self, state ):
        try:
            #print 'is', state, 'a bindings?'
            state( self._system._time_variable )
            #print 'yes!'
            # if it's a bindings, use it as is
        except TypeError:
            try:
                #print 'nope, is it a dict?'
                state[self._system._time_variable]
                # if it's a dict, make it a Bindings
                state = Bindings( state )
                #print 'yes, now it\'s', state
            except TypeError:
                #print 'nope, it must be a list'
                # otherwise, assume it's a list and make binding the default way
                state = Bindings( dict( (k,v) for k,v in
                    zip( [self._system._time_variable]+self._system._vars, state ) ) )
                #print 'now it\'s', state
        return state
    def append(self, state):
        """state: a dictionary to be added to self._timeseries"""
        #print 'append', state,
        state = self.point_to_bindings(state)
        #print 'as', state
        self._timeseries += [ self.point_to_bindings( state ) ]
        return self
    def __add__(self, other):
	# you can combine trajectories from different models, at your own risk
	#if self._system is not other._system:
	#    raise ValueError, "Can't concatenate incompatible trajectories"
	return Trajectory( self._system, self._timeseries + other._timeseries )
    def make_points(self, xexpr, yexpr):
        # for branching stuff
        # list the evaluations of the expressions at the points,
        # only at the ones where they both evaluate to numbers
        points = []
        xexpr = self._system._bindings(xexpr)
        yexpr = self._system._bindings(yexpr)
        for p in self._timeseries:
            try:
		print p,':',
                tup = (N(p(xexpr)), N(p(yexpr)))
		print tup
                points += [ tup ]
            except TypeError: pass
        return points
    def bind_value(self, ts, expr):
	try:
	    return ts(expr)
	except ValueError: # for instance, division by zero
	    return NaN
    def values(self, expr):
	bex = self._system._bindings(expr)
	return [ self.bind_value( t, bex ) for t in self._timeseries ]
    def plot(self, xexpr, yexpr, filename='', xlabel=-1, ylabel=-1, **args):
        """Make a 2-d plot of some pair of symbolic expressions that can
        be resolved to values of the state variables, time variable and
        parameters, along this trajectory.

        xexpr: expression to plot on the x axis.
        yexpr: expression to plot on the y axis.
        filename: filename to receive the plot.
        xlabel, ylabel: axis labels, if other than the text of the expressions.
        args: arguments to pass through to list_plot()
        """
        #print 'plot %s vs. %s' % (str(xexpr), str(yexpr))
        #print 'bindings are ', self._system._bindings
        xexpr = self._system._bindings(xexpr)
        try:
            if isinstance( yexpr, str ): raise TypeError
            #if isinstance( yexpr, basestr ): raise TypeError
            # TODO: this is not doing colors right
            P = None
	    colors = rainbow( len( yexpr ) )
            for y in yexpr:
                p = self.plot( xexpr, y, color=colors.pop(), legend_label='$%s$'%latex(y), **args )
                if P is None:
                    P = p
                else:
                    P += p
        except TypeError:
            yexpr = self._system._bindings(yexpr)
            #print 'after substitution: %s vs. %s' % (str(xexpr), str(yexpr))
            #print 'timeseries:', self._timeseries
            P = list_plot(
              self.make_points( xexpr, yexpr ),
              plotjoined = True,
              **args
            )
        if (xlabel == -1): xlabel = xexpr
        if (ylabel == -1): ylabel = yexpr
	try: basestring
	except NameError: basestring=str
	if not isinstance( xlabel, basestring ):
            xlabel = '$%s$' % latex(xlabel)
	if not isinstance( ylabel, basestring ):
            ylabel = '$%s$' % latex(ylabel)
        P.axes_labels( [xlabel,ylabel] )
        if (filename != ''):
            P.save(filename)
        return P
    def write_csv(self, filename, *columns):
        colnames = []
        for c in columns:
            try:
                c.expand() # see if it's an expression
            except AttributeError: # if not
                c = symbolic_expression(c) # it is now
            colnames = colnames + [c]
        import csv
        with open(filename, 'wb') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow( colnames )
            for dic in self._timeseries:
                csvwriter.writerow( [ c.substitute(dic) for c in colnames ] )

class Bindings(dict):
    """Definitions for variables and functions.

    When you have an ODE system, for example, du/dt = a - bu^2,
    and you apply the bindings {a:1, b:2}, you get du/dt = 1 - 2u^2.
    The former is better for working with symbolically, but you need
    the latter to integrate and get a numeric trajectory.

    We also represent equilibria using bindings, like {hat(u): sqrt(2)/2}.
    In the future we may use them for the ODE's flow as well, e.g. {u: 1-2*u^2}.
    """
    def __init__(self, *args, **named_args):
        #print 'Bindings:__init__', ('args = ', args,
        #    ', named_args =', dict(**named_args))
        self._function_bindings = FunctionBindings()
        for a in args:
            try: # works if it's a FunctionBindings
                a.merge_into_function_bindings(self)
            except AttributeError:
                try:
                    self.update(a) # works if it's a dict or a Bindings
                    try: # and in case it's a bindings
                        self._function_bindings.update( a._function_bindings )
                    except AttributeError: # ok, it wasn't a Bindings
                        pass
                except ValueError:
                    print 'Unrecognized initializer for bindings:', a
        self.update(dict(**named_args))
        for k,v in self.items():
            kk = symbolic_expression(k)
            if (kk != k):
                del self[k]
            self[kk] = symbolic_expression(v)
        #print ' =>', self
    def __repr__(self):
        frepr = repr(self._function_bindings)
        if frepr != '':
            frepr = ', ' + frepr
        return '{%s%s}' % (self.inner_repr(), frepr)
    def inner_repr(self):
        return ', '.join( '%s %s %s'%(k, '->', v) for k,v in self.items() ) 
    def _latex_(self):
        return '\\begin{align*}\n%s\n\\end{align*}' % self.latex_inner()
    def latex_inner(self):
        try:
            flx = self._function_bindings.latex_inner()
        except:
            flx = ''
        return '%s%s' % ( ' \\\\\n'.join( '  %s &\\to %s' % (latex(key), latex(val)) for key, val in self.items() ), flx )
    def substitute(self, expr):
        """Apply the bindings to an expression"""
	# if it has a bind method, use that
	try:
	    return expr.bind( self )
	except AttributeError: pass
        # if it's a Bindings or dict, apply ourself to all the values
        try:
            return expr.__class__( { k:self.substitute(v) for k,v in expr.items() } )
        except (AttributeError, ValueError): pass
        # or if it's a vector or something, apply to all entries
        try:
            return expr.apply_map( lambda x : self.substitute( x ) )
        except AttributeError: pass
        # in case expr is a plain string
        expr = symbolic_expression(expr)
        expr = expr.substitute_expression(self)
        #print ' => ', expr
        return self._function_bindings.substitute(expr)
    def __call__(self, expr):
        """Apply the bindings to an expression"""
        return self.substitute(expr)
    def __deepcopy__(self, _dict):
        """I'm not getting the results I want from deepcopy(bindings) --
        it seems to garble the function bindings' values -- so here's my own"""
        other = Bindings()
        other.update( self )
        other._function_bindings = deepcopy( self._function_bindings, _dict )
        return other
    # todo: merge() is akin to ODEsystem's bind(): it produces a bound copy
    # while bind() is akin to ODEsystem's bind_in_place(): it modifies self
    # this is problematic
    # in fact, unacceptable. fixed.
    def bind_in_place(self, *args, **named_args):
        print 'merge bindings:', self, ',', list(*args), ', ', dict(**named_args)
        other = Bindings(*args, **named_args)
        #for k, v in self.items():
        #    self[k] = other.substitute(v)
        try: # is it a FunctionBindings?
            other.merge_into_function_bindings(self)
        except AttributeError: # if not
            self.update(other)
            other._function_bindings.merge_into_function_bindings(self)
        self.apply_to_self()
        return self
    def bind(self, *args, **xrgs):
	return deepcopy(self).bind_in_place( *args, **xrgs )
    def merge(self, other={}, **args):
        """Combine with another set of bindings.  We assume that self is the
        bindings that have already been applied, and the other bindings are
        being applied afterward.  Thus self's bindings take priority, if there's
        any potential conflict."""
        return deepcopy(self).bind_in_place(other, **args)
    def __add__(self, other):
        return self.merge(other)
    def apply_to_self(self):
        """after merging bindings together, we have to apply the bindings to
        each other so that all substitutions get done in a single pass"""
        for k,v in self.items():
            self[k] = self.substitute(v)
        self._function_bindings.apply_bindings(self)

# see http://trac.sagemath.org/ticket/17553
# limit() and substitute_function() don't play well together.
# Once maxima returns a formal limit, it never gets re-evaluated
# even when the limit can be easily taken.  This re-evaluates them.
limop = limit( SR('f(x)'), x=0 ).operator()
def simplify_limits( expr ):
    return expr.substitute_function( limop, maxima_calculus.sr_limit )

from sage.symbolic.function_factory import function
class FunctionBindings(Bindings):
    """I am annoyed that we have to substitute functions in a different way
    from other things.  I hope to find a better way, by formalizing the idea
    that parameters sometimes depend on other things."""
    def __init__(self, *args, **named_args):
        for a in list( args ) + [ named_args ]:
	    for k,v in a.items():
                # don't store as before-and-after functions due to
                # function-pickling woes
                # http://trac.sagemath.org/ticket/17558
                # store as (function name, argument list): return value
                # i.e. with no function objects stored
                try:
                    # if it's a FunctionBindings or other dict { (<name>,<args>): <expr> }
                    # it's tricky to exclude string-valued ks from this
                    def rterr(): raise TypeError
                    self.update( { (isinstance(k,tuple) and k or rterr()):v } )
                except (TypeError, ValueError):
                    try:
                        # if it's a dict { <fn or name>: <fn or expr> }, including
                        # the kind you would use with substitute_function(),
                        # transform it to { (<name>,<args>):<expr> }
                        self.update( { (str(k),SR(v).arguments()):SR(v) } )
		    except TypeError:
		        # this comes up if it's a Piecewise object
		        # not sure what else
		        self.update( { (str(k),()):v } )
                    except ValueError:
                        print 'Unrecognized initializer for bindings:', {k:v}
    def __repr__(self):
        # would like to use unicode arrow u'\u2192' but causes output codec error
        return ', '.join( '%s(%s) %s %s' % (key[0], ','.join( str(k) for k in key[1] ), '->', str(val)) for key,val in self.items() )
    def latex_inner(self):
        return '\\\\\n'.join( '  %s(%s) &\\to %s' % ( latex(key[0]), ','.join( latex(a) for a in key[1] ), latex( val ) ) for key, val in self.items() )
    def substitute(self, expr):
        """Apply the bindings to an expression"""
        try:
            return expr.apply_map( lambda x : self.substitute( x ) )
        except AttributeError: pass
        #print 'substitute ', expr,' + ', self
        for k,v in self.items():
            # evade confusion between function('f', nargs=1)
            # and function('f'), since they get interchanged by
            # maxima and pickle operations
            # http://trac.sagemath.org/ticket/17503
	    try: expr = expr.substitute_function( function(k[0]), v.function( *k[1] ) )
	    except AttributeError: pass
            try: expr = expr.substitute_function( function(k[0], nargs=len(k[1])), v.function( *k[1] ) )
	    except AttributeError: pass
            expr = expr.substitute_function(function(k[0]),v)
            expr = expr.substitute_function(function(k[0],nargs=len(k[1])),v)
            #print ' => ', expr
        return simplify_limits( expr )
    def __deepcopy__(self, _dict):
        return FunctionBindings( self )
    def merge(self, other={}, **args):
        return Bindings(self).bind(other, **args)
    def merge_into_function_bindings(self, other):
        """Merge this FunctionBindings into other bindings, in place.
	
        self is a FunctionBindings while other is a regular Bindings.
        this function provides duck typing because Bindings doesn't have it"""
        #print 'merge',self,'into',other
        other._function_bindings.update(self)
    def apply_bindings(self, bindings):
        for k, v in self.items():
            self[k] = bindings(v)

class indexer(SageObject):
    """A little object that behaves like a dictionary of variables,
    but generates 'entries' according to a formula instead of storing
    them.  Used to fill in things like b_1 for population models.

    f: if f is a string, the indexer constructs an expression "f_i".
       if f is a function, the indexer uses the value of f(i).
       In either case, if the result is a string, it parses it into
       a sage expression object."""
    def __init__(self, f):
        self._f = f
    def __getitem__(self, i):
        try:
            return symbolic_expression(self._f(i))
        except TypeError:
            try:
                return self._f(i)
            except TypeError:
                return symbolic_expression("%s_%s" % (self._f,i))

def scriptedsymbol( base, superscripts=(), subscripts=() ):
    try:
        base.expand() # see if it's an expression
    except AttributeError: # if not
        base = SR.symbol( base ) # it is now
    name, latex_name = str(base), '{%s}'%latex(base)
    if len(superscripts) > 0:
        name += '^' + '^'.join(str(s) for s in superscripts)
        if len(superscripts) > 1:
            latex_name += '^{%s}' % ''.join('{%s}'%latex(s) for s in superscripts)
        else:
            latex_name += '^{%s}' % latex(superscripts[0])
    if len(subscripts) > 0:
        name += '_' + '_'.join(str(s) for s in subscripts)
        if len(subscripts) > 1:
            latex_name += '_{%s}' % ''.join('{%s}'%latex(s) for s in subscripts)
        else:
            latex_name += '_{%s}' % latex(subscripts[0])
    return SR.symbol( name, latex_name=latex_name )

def subscriptedsymbol( base, *subscripts ):
    return scriptedsymbol( base, subscripts=subscripts )

# this class was nested inside indexer_2d, but it broke pickling
class indexer_2d_inner(indexer):
    def __init__(self, f, i):
	self._f = f
	self._i = i
    def __getitem__(self, j):
	return SR.symbol( '%s_%s_%s' % (self._f, self._i, j),
	    latex_name='{%s}_{%s%s}' % (self._f, self._i, j) )

class indexer_2d(indexer):
    """Instead of mapping i |-> x_i, this does a 2-step mapping
    i |-> j |-> x_i_j.  That is, indexer_2d('x')[i][j] produces x_i_j."""
    def __getitem__(self, i):
        return indexer_2d_inner( self._f, i )

class indexer_2d_reverse_inner(indexer):
    def __init__(self, f, j):
	self._f = f
	self._j = j
    def __getitem__(self, i):
	return SR.symbol( '%s_%s_%s' % (self._f, i, self._j ),
	    latex_name='%s_{%s%s}' % (self._f, i, self._j ) )

class indexer_2d_reverse(indexer_2d):
    """Just like indexer_2d but maps j |-> i |-> x_i_j rather than to
    x_j_i as the other one would do.  Useful when you want a way to
    generate all x_i_j for a given j, as I do."""
    def __getitem__(self, j):
        return indexer_2d_reverse_inner( self._f, j )

class const_indexer(indexer):
    """always return the same value"""
    def __getitem__(self, i):
        return self._f

# And now the dynamical systems classes.

class ODESystem(SageObject):
    """A system of ordinary differential equations.

    This object holds a symbolic representation of the ODE's flow vector field,
    and can use it in numerical integration as well as making it available
    for symbolic manipulations.
    """
    def __init__(self, flow, vars, time_variable=SR.var('t'),
            bindings=Bindings()):
        """Construct from the basic data.

        flow: a dictionary with one entry for each state variable.
            Keys are variables, and values are symbolic expressions.
        vars: a list of the state variables.  This is conceptually
            equivalent to flow.keys(), but it's important to give them
            an ordering.  These are variables, not strings.
        time_variable: typically t, but could be something else.
        bindings: optional substitutions of values for symbolic variables
            and functions.  This is used to interpret symbolic
            expressions for plotting, etc.
        """
        #print 'ODESystem init', (flow, vars, time_variable, bindings)
        self._flow = flow
        self._vars = vars
        self._time_variable = time_variable
        self._bindings = bindings
        self._flow = { k:bindings(v) for k,v in self._flow.items() }
    def __repr__(self):
        """Output the system as a system of differential equations"""
        return join( ('%s -> %s'%(v,self._flow[v])
            for v in self._vars), '\n')
    def x__deepcopy__(self, _dict):
        """There seems to be a weird bug when you deepcopy a symbolic expression.
        So I'm trying to avoid that."""
        other = copy( self )
        other._flow = dict()
        other._flow.update( self._flow )
        other._bindings = deepcopy( self._bindings )
        return other
    def _latex_(self):
        """Output the system as a system of differential equations in LaTeX form"""
        return self.latex_text() # not correct in math mode!
    def latex_text(self):
        return ('\\iflatexml\n' +
	    '\\begin{align*}\n' + '\\\\\n'.join(
            r'\frac{d%s}{d%s} &= %s'%(latex(v),latex(self._time_variable),latex(self._flow[v]))
                for v in self._vars ) + '\n\\end{align*}\n'
            '\\else\n' +
	    '\\begin{dgroup*}\n\\begin{dmath*}\n' + '\\end{dmath*}\n\\begin{dmath*}\n'.join(
             r'\frac{d%s}{d%s} = %s'%(latex(v),latex(self._time_variable),latex(self._flow[v]))
                for v in self._vars ) + '\n\\end{dmath*}\n\\end{dgroup*}\n' +
	    '\\fi\n')
    def write_latex(self, filename, inline=False):
        """Output the system in LaTeX form to a file"""
	if inline:
            ltxout = open( filename, 'w' )
	else:
	    import latex_output
            ltxout = latex_output.latex_output( filename )
        ltxout.write( self.latex_text() )
        ltxout.write( '\n\\vspace{24pt}\n' )
        ltxout.close()
    def flow( self, at=None ):
	if at is None: at = self._vars
	try:
	    at( self._vars[0] )
	except TypeError: # if it's not a Bindings, make it one
	    at = Bindings( zip( self._vars, at ) )
	return vector( [ at(self._flow[v]) for v in self._vars ] )
    def time_variable(self):
        """Provide the independent variable, for instance for plotting"""
        return self._time_variable
    def bind_in_place(self, *bindings, **args):
        """Apply bindings to my formulas without copying to a new object.

        See bind(), below."""
        binding = Bindings( *bindings, **args )
        self._flow = { k:binding(v) for k,v in self._flow.items() }
        self._bindings = self._bindings + binding
    def bind(self, *bindings, **args):
        """If you create a system with various symbolic parameters, like
        dy/dt = ay^2 + by + c, or something, you can't numerically
        integrate it as is, you have to give values to parameters a, b, and c.
        This method gives you a system just like self, but with parameters
        bound to the values provided.

        bindings: Bindings objects recording values for some variables and/or
        functions.
        """
        bound = deepcopy( self )
        bound.bind_in_place( *bindings, **args )
        return bound
    def solve(self, initial_conditions, end_time=20, start_time=0, step=0.1):
	return self.desolve( initial_conditions, start_time=start_time, end_time=end_time, step=step)
    def desolve(self, initial_conditions, end_time=20, start_time=0, step=0.1):
        """Use a numerical solver to find a concrete trajectory of the system.

        initial_conditions: list or Bindings of initial values for the state variables"""
	## if initial_conditions is a Bindings, make a list
	try: initial_conditions = [ initial_conditions(x) for x in self._vars ]
	except TypeError: 
	    ## or if it's a dict, make a list
	    try: initial_conditions = [ initial_conditions[x] for x in self._vars ]
	    except TypeError: pass
        print "desolve: %s, %s, %s, ivar=%s, end_points=%s, step=%s" % (
          [self._flow[v] for v in self._vars],
          self._vars,
          [start_time] + initial_conditions,
          self._time_variable,
          end_time,
          step )
	print "initial flow:", [ self._flow[x].subs( **( { str(k):v for k,v in zip([self._time_variable] + self._vars, [start_time] + initial_conditions ) } ) ) for x in self._vars ]
        soln = desolve_system_rk4(
          [self._flow[v] for v in self._vars],
          self._vars, [start_time] + initial_conditions,
          ivar=self._time_variable,
          end_points=end_time, step=step )
        #return Trajectory(self, soln)
	# this seems to help when sympy is involved
	return Trajectory(self, [ [N(z) for z in l] for l in soln ])
    def plot_vector_field(self, xlims, ylims, filename='', vf=None, xlabel=-1, ylabel=-1, **args):
        xlims = tuple( self._bindings.substitute( v ) for v in xlims )
        ylims = tuple( self._bindings.substitute( v ) for v in ylims )
        xvar, yvar = ( xlims[0], ylims[0] )
        if vf is None:
            vf = [ self._flow[v] for v in [ xvar, yvar ] ]
        PL = plot_vector_field( vf, xlims, ylims, **args )
        if (xlabel == -1): xlabel = xvar
        if (ylabel == -1): ylabel = yvar
        PL.axes_labels( ['$%s$'%latex(v) for v in (xlabel,ylabel)] )
        if (filename != ''):
            PL.save(filename)
        return PL
    def plot_isoclines(self, xlims, ylims, isoclines=None, filename='', xlabel=-1, ylabel=-1, **args):
        xlims = tuple( self._bindings.substitute( v ) for v in xlims )
        ylims = tuple( self._bindings.substitute( v ) for v in ylims )
        if isoclines is None:
            isoclines = [ (v,0) for v in self._vars ]
        p = Graphics()
        for v, z in isoclines:
            p += implicit_plot( self._flow[v] == z, xlims, ylims, **args )
        if (xlabel == -1): xlabel = xlims[0]
        if (ylabel == -1): ylabel = ylims[0]
        p.axes_labels( ['$%s$'%latex(v) for v in (xlabel,ylabel)] )
        if (filename != ''):
            p.save(filename)
        return p
    def add_hats(self, expression=None):
        """convert e.g. X_i to Xhat_i.  Should only
        operate on the state variables."""
        try:
            self._add_hats( self._time_variable )
        except (AttributeError,TypeError):
            self._add_hats = Bindings( { v:hat(v) for v in self._vars } )
        # self.add_hats( expression ) returns hat-added expression
        if expression is not None:
            return self._add_hats( expression )
        # self.add_hats() returns the add-hats bindings
        return self._add_hats
    def remove_hats(self):
        return Bindings( { v:k for k,v in self.add_hats().items() } )
    def equilibrium_vars(self):
        add_hats = self.add_hats()
        return [ add_hats(v) for v in self._vars ]
    def symbolic_equilibria( self ):
	return self.add_hats()
    def limit( self, **lims ):
	"""limit( var=val, ... ):
	Return a transformed dynamical system, equal to this system with the
	specified limits applied.
	Can perform multiple transformations including
	* Separate fast and slow timescales by changing parameters to
	  zero, infinite values, or other threshold values
	* Extract subsystems or other limiting cases by changing state
	  variables to zero, or infinite or other values.
	These functionalities will be implemented as needed."""
	return deepcopy(self).limit_in_place( **lims )
    def limit_in_place( self, **lims ):
	for k,v in lims.iteritems():
	    ksr = SR(k)
	    if ksr.is_symbol():
		self._vars = [ x for x in self._vars if x != ksr ]
		## to do: specialized version of this in BoxModel, to
		## make transitions collapse to identity of boxes, and
		## to change finite boxes to infinite sources, etc.
		self._flow = { x:self._flow[x].limit( **{k:v} ) for x in self._vars }
		self._bindings[ksr] = v
	    else:
		raise ValueError, "Not skillful enough to take limit {0} -> {1}".format( str(k), str(v) )
	return self
    def equilibria(self, *ranges, **opts):
	solve_numerically = opts.get('solve_numerically',False)
        add_hats = self.add_hats()
	equil_vars = self.equilibrium_vars()
	if solve_numerically:
	    ranges = { l[0]:l[1:] for l in ranges }
            def grid_range( var ):
                try:
                    return numpy.arange( *ranges[var] )
                except (TypeError, AttributeError):
                    return numpy.arange( -1, 1.01, 0.5 )
            grid_ranges = ( grid_range(var) for var in self._vars )
	    flows = [ add_hats(rhs) for rhs in self._flow.values() ]
	    def flowfun(vec):
		dic = { x:v for x,v in zip(equil_vars,vec) }
		return [ f.subs( dic ) for f in flows ]
	    sols = set()
	    import itertools, scipy.optimize
            for grid_point in itertools.product( *grid_ranges ):
		print grid_point; sys.stdout.flush()
		res = scipy.optimize.root( flowfun, grid_point, tol=1e-10 )
		## caution: res.success is reporting false positives!
		if res.success and all( round(f,ndigits=5)==0 for f in res.fun ):
	            sols.add( tuple( [ round(N(x),ndigits=5) for x in res.x ] ) )
	    equilibria = [ dict( zip( equil_vars,sol ) ) for sol in sols ]
	else:
            equil_eqns = [ 0 == add_hats(rhs) for rhs in self._flow.values() ]
            equilibria = solve( equil_eqns, *self.equilibrium_vars(), solution_dict=True )
        #equilibria = [ { k:(v) for k,v in soln.items() } for soln in solns ]
        return equilibria
    def interior_equilibria(self):
        return [ eq for eq in self.equilibria() if all( eq[hat(x)] != 0 for x in self._vars ) ]
    def jacobian(self, at=None):
        j = dict( ( ki, dict( ( kj, diff( self._flow[ki], kj ) ) for kj in self._vars ) ) for ki in self._vars )
        try:
            j = dict( ( ki, dict( ( kj, v.substitute( at ) ) for (kj, v) in d.items() ) ) for (ki, d) in j.items() )
        except:
            pass
        return j
    def jacobian_matrix(self, at=None):
        jdict = self.jacobian(at)
        return matrix( [ [ jdict[ki][kj] for kj in self._vars ] for ki in self._vars ] )
    def is_stable(self, matrix):
        """In an ODE system stability means Re(eigenvalue) < 0.  In a
        difference equation it would be different..."""
        for e in matrix.eigenvalues():
            if real( e ) >= 0: return false
        return true
    def stable_equilibria(self):
        remove_hats = self.remove_hats()
        return [ eq for eq in self.equilibria() if self.is_stable( self.jacobian_matrix( { remove_hats(k):v for k,v in eq.items() } ) ) ]
    def stable_interior_equilibria(self):
        remove_hats = self.remove_hats()
        return [ eq for eq in self.interior_equilibria() if self.is_stable( self.jacobian_matrix( { remove_hats(k):v for k,v in eq.items() } ) ) ]
    def test_equilibrium( self, eq, t=0 ):
        '''evaluate whether a state vector is an equilibrium (or very near).

        eq: an array, vector, tuple, or Bindings.
        returns True or False'''
        try: # is it a Bindings?
            eq('t')
        except AttributeError: # if eq() doesn't work - is this right?
            eq = Bindings( zip(self._vars, eq) )
        flow = self._flow( eq(v) for v in self._vars )
        return ( sum( x*x for x in flow ) < 0.001 )
    def plot_bifurcation_diagram( self, ind_range, dep_range, filename=None, xlabel=None, ylabel=None, **args ):
        (p, pmin, pmax) = ind_range
        (y, ymin, ymax) = dep_range
        def bifurcation_points( fixed_points_eigenvalues, p ):
            '''Given a list of eigenvalues, and range of "p" values,
            find values of parameter "p" where dominant eigenvalue changes sign
            fixed_points_eigenvalues: list of lists: a list of eigenvalues for each fixed point
            p: the symbolic variable to vary as the independent parameter'''
            # todo: teach population dynamics class about physically admissible equilibria?
            bif_points = set()
            for ev in fixed_points_eigenvalues:
                # ev is the eigenvalues at one fixed point
                # find zeros of all those
                possible_p = [ solve( v == 0, p, solution_dict=true ) for v in ev ]
                # and also their poles
                possible_p += [ solve( 1/v == 0, p, solution_dict=true ) for v in ev ]
                # that unfortunately is a list of lists of dictionaries
                possible_p = [ ps for pp in possible_p for ps in pp ]
                # now it's a list of dictionaries
                possible_p = [ pp[p] for pp in possible_p ]
                # now it's a list of the actual values of the parameter
                # those are places where any eigenvalue might change sign
                bif_points = bif_points.union( bp for bp in set( possible_p )
                    if ( max( v.substitute( { p : bp - 0.0001 } ) for v in ev ) *
                        max( v.substitute( { p : bp + 0.0001 } ) for v in ev ) < 0 ) )
                # now these are places where the dominant eigenvalue changes sign
            return list( bif_points )
        equilibria_b = self.equilibria()
        remove_hats = self.remove_hats()
        jacobians_b = [ self.jacobian_matrix( { remove_hats(k):v for k,v in eq.items() } ) for eq in equilibria_b ]
        eigenvalues_b = [ j.eigenvalues() for j in jacobians_b ]
        bif_points = bifurcation_points( eigenvalues_b, p )
        # plot!
        bp = Graphics()
        for s, e in zip( [pmin] + bif_points, bif_points + [pmax] ):
            for eq, ev in zip( equilibria_b, eigenvalues_b ):
                vmid = [ v.substitute( { p : (s + e)/2 } ) for v in ev ]
                dom_ev = max( vmid )
                bp += plot( y.substitute( eq ), (p, s, e), color =
                    ('black' if dom_ev < 0 else 'gray' if dom_ev == 0
                       else 'blue' if min( vmid ) < 0 else 'red'),
                    thickness=2, ymin=ymin, ymax=ymax, **args )
        if xlabel is None: xlabel = '$%s$' % latex( p )
        if ylabel is None: ylabel = '$%s$' % latex( y )
        bp.axes_labels( [ xlabel, ylabel ] )
        if filename is not None:
            bp.save( filename )
        return bp
    def gsl_system(self):
        il = range(len(self._vars))
        cython_code = (
            "cimport sage.gsl.ode\n"
            "import sage.gsl.ode\n"
            "from sage.ext.interpreters.wrapper_rdf cimport Wrapper_rdf\n"
            "include 'gsl.pxi'\n"
            "cdef class gsl_ode_system(sage.gsl.ode.ode_system):\n"
            + ''.join(
            "    cdef Wrapper_rdf _flow%d\n"%i for i in il
            ) +
            "    def __init__( self, "
            + ', '.join("Wrapper_rdf f%d"%i for i in il)
            + " ):\n"
            + ''.join(
            "        self._flow%d = f%d\n"% (i,i) for i in il
            ) +
            "    cdef int c_f( self, double t, double *y, double *dydt ):\n"
            #"        cdef double *result = [t]\n"
            + (''.join(
            "        self._flow%d.call_c( y, &dydt[%d] )\n"
            #"        dydt[%d] = result[0]\n"
                % (i,i) for i in il
            )) +
            "        return GSL_SUCCESS\n"
            # need c_j() as well?
        )
        #print cython_code
        module = sage.misc.cython.compile_and_load( cython_code )
        T = ode_solver()
        #T.algorithm = 'bsimp'
        T.function = module.gsl_ode_system( *(
            fast_float( self._flow[v], *self._vars )
            for v in self._vars
        ) )
        return T
    def solve_gsl( self, initial_conditions, start_time=0, end_time=20, step=0.1):
        try:
            self._gsl_system
        except AttributeError:
            self._gsl_system = self.gsl_system()
        self._gsl_system.ode_solve(
            t_span=[start_time, end_time],
            num_points = (end_time - start_time) / step + 1,
            y_0 = initial_conditions
        )
        return Trajectory( self, ( [ t ] + ys for t,ys in self._gsl_system.solution ) )
        return self._gsl_system.solution

# used by NumericalODESystem.solve()
# unlike the standard odeint(), this fails when the du/dt function
# raises and exception
# http://www.pythoneye.com/212_16439560/
import scipy.integrate, numpy
def fake_odeint(system, y0, t, args=None, Dfun=None):
    def cf(t, y, s):
        scf = list( system.compute_flow(y, t, True) )
        return scf
    ig = scipy.integrate.ode(cf, Dfun)
    ig.set_integrator('lsoda', method='adams')
    ig.set_initial_value(y0, t=t[0])
    if not isinstance(args, (tuple,list)): args = (args,)
    ig.set_f_params(*args)
    y = [y0]
    for tt in t[1:]:
        #if not ig.successful():
        #    break
        try:
            y.append(ig.integrate(tt))
        except (EquilibriumDetectedException,UnboundedDynamicsException) as ex:
            raise TrajectoryInterruptedException(
                Trajectory(
                    system,
                    [ [tt]+list(point) for tt,point in zip(t,y) ]
                ),
                ex
            )
    return numpy.array(y)

class NumericalODESystem( ODESystem ):
    """Where the base ODESystem has a system of Sage expressions to
    specify its flow, this uses an instance method, operating on
    python data (compatible with scipy.integrate.odeint) to provide
    the flow vector field.
    
    This is an abstract base class. It's used by providing a subclass
    with a meaningful compute_flow() method."""
    def __init__(self, vars, time_variable=SR.var('t'),
            bindings=Bindings(), flow=None):
        self._vars = vars
        self._time_variable = time_variable
        self._bindings = bindings
        # note flow is only used for printing out the system
        self._flow = flow
        if self._flow is not None:
            self._flow = { k:bindings(v) for k,v in self._flow.items() }
    def solve(self, initial_conditions, start_time=0, end_time=20, step=0.1):
	print 'NumericalODESystem.solve'
        times = numpy.arange( start_time, end_time, step )
        soln = fake_odeint( self, numpy.array( initial_conditions, float ), times, self )
        # make that timeseries of lists into timeseries of dictionaries
        timeseries = [ dict( (k,v) for k,v in
          zip([self._time_variable]+self._vars,[t]+list(point)) )
          for t, point in zip(times, soln) ]
        return Trajectory(self, timeseries)
    def compute_flow(self, y, t):
        # provide this in a subclass.  It needs to return a numpy array of
        # floats providing the dy/dt vector matching the y argument.
        pass
    def detect_equilibrium( self, eq, flow=None, t=0 ):
        '''evaluate whether a state vector is an equilibrium (or very near).

        eq: an array, vector, tuple, or Bindings.
        returns True or False'''
        try: # in case it's a Bindings
            eq = numpy.array( [ eq(v) for v in self._vars ] )
        except TypeError: # if eq() doesn't work - is this right?
            pass
        if flow is None:
            flow = self.compute_flow( eq, t )
        flow_norm = sum( x*x for x in flow )
        #print 'check', flow_norm, 'vs 0.0005'
        # are we close to equilibrium
        if flow_norm < 0.0001:
            # yes
            #print 'norm is small:',flow_norm
            # but are we getting closer
            next_state = numpy.array( [ x + 0.001*dx for x, dx in zip( eq, flow ) ], float )
            next_flow = self.compute_flow( next_state, t )
            next_norm = sum( x*x for x in next_flow )
            #print 'next norm is', next_norm
            if next_norm < flow_norm:
                # yes
                #print 'equilibrium!'
                return True
            # no
        return False
    def equilibria(self, ranges=None):
        def grid_range( var ):
            try:
                return ranges[var]
            except (TypeError, AttributeError):
                return numpy.arange( -1, 1.01, 0.5 )
        grid_ranges = ( grid_range(var) for var in self._vars )
        import sets, itertools
        eq_set = sets.Set()
        def objective_function( x ):
            #print 'try at', x, '...',
            try:
                answer = sum( d*d for d in self.compute_flow( x, 0 ) )
                #print answer
                return answer
            except (ValueError, DynamicsException), e:
                #print 'nope'
                return oo
        for grid_point in itertools.product( *grid_ranges ):
            #print 'grid point', grid_point
            sys.stdout.flush()
            try: # check the flow exists at the start point
                self.compute_flow( grid_point, 0 )
            except (DynamicsException,ValueError):
                continue
            #print 'minimize'
            sys.stdout.flush()
            minimum = minimize( objective_function, grid_point )
            if objective_function( minimum ) < 0.001:
                minimum = tuple( N(m, 10) for m in minimum )
                #print 'found', minimum
                sys.stdout.flush()
                eq_set.add( minimum )
        add_hats = self.add_hats()
        equilibria = [ Bindings( { add_hats(k):v for k,v in zip( self._vars, zero ) } ) for zero in eq_set ]
        return equilibria
    def plot_vector_field(self, xlims, ylims, filename='', vf=None, xlabel=-1, ylabel=-1, t=0, **options):
        from sage.plot.all import Graphics
        from sage.misc.misc import xsrange
        from sage.plot.plot_field import PlotField
        from sage.plot.misc import setup_for_eval_on_grid
        xlims = tuple( self._bindings.substitute( v ) for v in xlims )
        ylims = tuple( self._bindings.substitute( v ) for v in ylims )
        xvar, yvar = ( SR(xlims[0]), SR(ylims[0]) )
        if vf is None:
            vf = [ self._flow[v] for v in [ xvar, yvar ] ]
        _, rangez = setup_for_eval_on_grid( [ xvar, yvar ], ( xlims, ylims ), options['plot_points'] )
        xs, ys, dxs, dys = [], [], [], []
        for x in xsrange( *rangez[0], include_endpoint=True ):
            for y in xsrange( *rangez[1], include_endpoint=True ):
                try:
                    v = self.compute_flow( numpy.array( [ x, y ], float ), t )
                    xs.append( x )
                    ys.append( y )
                    if -0.001 <= v[0] <= 0.001: v[0] = 0
                    dxs.append( v[0] )
                    if -0.001 <= v[1] <= 0.001: v[1] = 0
                    dys.append( v[1] )
                except: pass
        dxs = numpy.ma.masked_invalid( numpy.array( dxs, dtype=float ) )
        dys = numpy.ma.masked_invalid( numpy.array( dys, dtype=float ) )
        g = Graphics()
        g._set_extra_kwds( Graphics._extract_kwds_for_show( options ) )
        if len(xs) > 0:
            g.add_primitive( PlotField( xs, ys, dxs, dys, options ) )
        if (xlabel == -1): xlabel = xvar
        if (ylabel == -1): ylabel = yvar
        g.axes_labels( ['$%s$'%latex(v) for v in (xlabel,ylabel)] )
        return g

class PopulationDynamicsSystem(ODESystem):
    """Specialization of ODESystem for population dynamics systems"""
    def __init__(self, vars, population_indices, population_indexer,
            time_variable=SR.var('t'),
            bindings=Bindings()):
        """for the purposes of AdaptiveDynamicsModel, we need our population
        dynamics models to be fairly general, and able to produce a specific
        flow for a varying number of "population indices".

        population_indices is a list of indices, for instance [1,2,3].
        population_indexer makes state variables from those indices, so, for
        example, population_indexer[0] should be a variable, say, X_0.
        flow() should generate the model's dynamics given the current
        population_indices and whatever else it needs.  It should be able to
        be called more than once."""
        self._nonpop_vars = vars
        super(PopulationDynamicsSystem, self).__init__({},
            vars, time_variable, bindings)
        self._population_indexer = population_indexer
        # do this last because it computes the flow.
        self.set_population_indices(population_indices)
    def population_vars(self):
        return [self._population_indexer[i] for i in self._population_indices]
    def n_populations(self):
	return len( self._population_indices )
    # don't set _population_indices directly, call this, to keep the flow in sync
    def set_population_indices(self, xi):
	print 'set_population_indices:', xi
        self._population_indices = xi
        self._vars = self._nonpop_vars + self.population_vars()
        self._flow = self.flow()
        #self._flow = dict( (k,self._bindings(v)) for k,v in self._flow.items() )
        self._add_hats = None
    def remove_population(self, i):
	self.set_population_indices( [ j for j in self._population_indices if j != i ] )
    def remove_populations(self, ps):
	sps = set(ps)
	self.set_population_indices( [ j for j in self._population_indices if j not in sps ] )
    # subclasses have to provide the flow
    def flow(self): pass
    def mutate(self, resident_index):
        mutant_index = 1 + max( self._population_indices )
        self.set_population_indices( self._population_indices + [ mutant_index ] )
        return mutant_index
    def fake_population_index(self):
        return 'x'
    def nontrivial_equilibria(self):
        equilibria = self.equilibria()
        return [ eq for eq in equilibria if sum( eq[hat(x)] for x in self.population_vars() ) != 0 ]
    def stable_nontrivial_equilibria(self):
        remove_hats = self.remove_hats()
        return [ eq for eq in self.nontrivial_equilibria() if self.is_stable( self.jacobian_matrix( { remove_hats(k):v for k,v in eq.items() } ) ) ]
    def interior_equilibria(self):
        equilibria = self.equilibria()
        def is_interior(x):
            try: return N(x) > 0
            except TypeError: return x != 0
        return [ eq for eq in equilibria if all( is_interior(eq[hat(x)]) for x in self.population_vars() ) ]
    def stable_interior_equilibria(self):
        remove_hats = self.remove_hats()
        return [ eq for eq in self.interior_equilibria() if self.is_stable( self.jacobian_matrix( { remove_hats(k):v for k,v in eq.items() } ) ) ]
    def plot_ZNGIs( self, xlims, ylims, vars=None, **args ):
        """Plot populations' zero net growth isoclines, using ODESystem's
        plot_isoclines()."""
        if vars is None:
            vars = self.population_vars()
        return self.plot_isoclines( xlims, ylims, [ (v,0) for v in vars ], **args )