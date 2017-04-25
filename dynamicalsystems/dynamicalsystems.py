"""Dynamical systems classes"""
## these 2 lines seem to be needed for doctesting
## http://webcache.googleusercontent.com/search?q=cache:E3ULGhsSTZMJ:ask.sagemath.org/question/7998/trouble-using-import-and-doctest-together/+&cd=6&hl=en&ct=clnk&gl=us
import sys
sys.path.append(".")
##

from sage.all import *
from string import join
from sage.misc.latex import _latex_file_
latex.add_to_preamble('\\usepackage{amsmath}')

from latex_output import *
from bindings import *

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
    def make_points(self, xexpr, yexpr, zexpr=None):
        # for branching stuff
        # list the evaluations of the expressions at the points,
        # only at the ones where they both evaluate to numbers
        points = []
        xexpr = self._system._bindings(xexpr)
        yexpr = self._system._bindings(yexpr)
        zexpr = self._system._bindings(zexpr)
        for p in self._timeseries:
            try:
                #print p,':',
                if zexpr is not None:
                    tup = (N(p(xexpr)), N(p(yexpr)), N(p(zexpr)))
                else:
                    tup = (N(p(xexpr)), N(p(yexpr)))
                #print tup
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
    def plot(self, xexpr, yexpr, zexpr=None, filename='', xlabel=-1, ylabel=-1, label_transformation=lambda x:x, **args):
        """Make a 2-d or 3-d plot of some pair of symbolic expressions that can
        be resolved to values of the state variables, time variable and
        parameters, along this trajectory.

        xexpr: expression to plot on the x axis.
        yexpr: expression to plot on the y axis.
        zexpr: expression, if any, to plot on the z axis.
        filename: filename to receive the plot.
        xlabel, ylabel: axis labels, if other than the text of the expressions.
        (labels not implemented in 3-d.)
        args: arguments to pass through to list_plot()
        """
        #print 'plot %s vs. %s' % (str(xexpr), str(yexpr))
        #print 'bindings are ', self._system._bindings
        xexpr = self._system._bindings(xexpr)
        try:
            # catch the "y is a string" case, because the test below won't
            if isinstance( yexpr, str ): raise TypeError('yexpr is a string, "{}"'.format(yexpr))
            #if isinstance( yexpr, basestr ): raise TypeError
            # TODO: this is not doing colors right
            P = None
            # this or the next line raise TypeError if y is not a list or tuple
            colors = rainbow( len( yexpr ) )
            for y in yexpr:
                #print 'iterate, plot',y
                zargs = copy(args)
                if 'color' not in zargs: zargs['color'] = colors.pop() 
                if 'legend_label' not in zargs: zargs['legend_label'] = '$%s$'%latex( label_transformation(y) )
                py = self.plot( xexpr, y, **zargs )
                if P is None:
                    P = py
                else:
                    P += py
            ylabel = ''
        except TypeError, te:
            #print 'TypeError:', te
            yexpr = self._system._bindings(yexpr)
            #print 'after substitution: %s vs. %s' % (str(xexpr), str(yexpr))
            #sys.stdout.flush()
            #print 'timeseries:', self._timeseries
            if zexpr is not None:
                import sage.plot.plot3d.shapes2
                P = sage.plot.plot3d.shapes2.Line(
                  self.make_points( xexpr, yexpr, zexpr ),
                  **args
                )
            else:
                P = list_plot(
                  self.make_points( xexpr, yexpr ),
                  plotjoined = True,
                  **args
                )
        if zexpr is None:
            if (xlabel == -1): xlabel = xexpr
            if (ylabel == -1): ylabel = yexpr
            try: basestring
            except NameError: basestring=str
            if not isinstance( xlabel, basestring ):
                xlabel = '$%s$' % latex( label_transformation(xlabel) )
            if not isinstance( ylabel, basestring ):
                ylabel = '$%s$' % latex( label_transformation(ylabel) )
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
    def __inner_class(self):
        return indexer_2d_inner
    def __init__(self, f, inner_class=None):
        if inner_class is None: self._inner_class = self.__inner_class()
        else: self._inner_class = inner_class
        super(indexer_2d,self).__init__(f)
    def __getitem__(self, i):
        return self._inner_class( self._f, i )

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
    #def __getitem__(self, j):
    #    return indexer_2d_reverse_inner( self._f, j )
    def inner_class(self):
        return indexer_2d_reverse_inner

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
    def bind(self, *bindings, **args):
        """If you create a system with various symbolic parameters, like
        dy/dt = ay^2 + by + c, or something, you can't numerically
        integrate it as is, you have to give values to parameters a, b, and c.
        This method gives you a system just like self, but with parameters
        bound to the values provided.

        bindings: Bindings objects recording values for some variables and/or
        functions.
        args: key-value pairs pairing names to values
        """
        return deepcopy( self ).bind_in_place( *bindings, **args )
    def bind_in_place(self, *bindings, **args):
        """Apply bindings to my formulas without copying to a new object.
        See bind(), above.
        """
        b = Bindings( *bindings, **args )
        self._bindings.merge_in_place( b )
        self._flow = { k:self._bindings(v) for k,v in self._flow.items() }
        return self
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
        #print "desolve: %s, %s, %s, ivar=%s, end_points=%s, step=%s" % (
        #  [self._flow[v] for v in self._vars],
        #  self._vars,
        #  [start_time] + initial_conditions,
        #  self._time_variable,
        #  end_time,
        #  step )
        #print "initial flow:", [ self._flow[x].subs( **( { str(k):v for k,v in zip([self._time_variable] + self._vars, [start_time] + initial_conditions ) } ) ) for x in self._vars ]
        soln = desolve_system_rk4(
          [self._flow[v] for v in self._vars],
          self._vars, [ SR(x0) for x0 in [start_time] + initial_conditions ],
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
        return Bindings( { v:k for k,v in self.add_hats()._dict.items() } )
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
        substitutions = Bindings( opts.get('substitutions', {}) )
        add_hats = self.add_hats()
        equil_vars = self.equilibrium_vars()
        if solve_numerically:
            ranges = { l[0]:l[1:] for l in ranges }
            def grid_range( var ):
                try:
                    return numpy.arange( *ranges[var] )
                except (KeyError, TypeError, AttributeError):
                    return numpy.arange( -1, 1.01, 0.5 )
            grid_ranges = ( grid_range(var) for var in self._vars )
            flows = [ add_hats( substitutions( self._flow[v] ) ) for v in self._vars ]
            def flowfun(vec):
                dic = { x:v for x,v in zip(equil_vars,vec) }
                fvec = [ f.subs( dic ) for f in flows ]
                #print fvec
                return [ N( f ) for f in fvec ]
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
            equil_eqns = [ 0 == add_hats( substitutions( rhs ) ) for rhs in self._flow.values() ]
            equilibria = solve( equil_eqns, *self.equilibrium_vars(), solution_dict=True )
            #equilibria = [ eq for eq in equilibria if all( add_hats(f).subs(eq) == 0 for f in self._flow.values() ) ]
        return equilibria
    def nontrivial_equilibria(self):
        equilibria = self.equilibria()
        print equilibria
        try: return [ eq for eq in equilibria if any( eq[hat(x)] != 0 for x in self._vars ) ]
        except KeyError: return []
    def interior_equilibria(self):
        try: return [ eq for eq in self.equilibria() if all( eq[hat(x)] != 0 for x in self._vars ) ]
        except KeyError: return []
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
    def gsl_system(self, make_t_explicit=False, version=2):
	## compatibility with multiple Sage releases
	##  version == 0: pre 7.5?
	##  version == 1: 7.5? 7.6?
	##  version == 2: at least some 7.6 releases
        if make_t_explicit:
            il = range(len(self._vars)+1)
            flow = copy(self._flow)
            flow[ self._time_variable ] = SR(1)
            vars = self._vars + [ self._time_variable ]
            print 'making gsl system with t'; sys.stdout.flush()
        else:
            il = range(len(self._vars))
            flow = self._flow
            vars = self._vars
            print 'making gsl system without t'; sys.stdout.flush()
        cython_code = (
              (
	    "cimport sage.gsl.ode\n"
            "import sage.gsl.ode\n"
		if version < 2 else
	    "cimport sage.calculus.ode\n"
	    "import sage.calculus.ode\n"
	      ) +
            "from sage.ext.interpreters.wrapper_rdf cimport Wrapper_rdf\n"
	      + (
	    "include 'gsl/gsl.pxi'\n"
	        if version==0 else
            "from sage.libs.gsl.all cimport *\n"
	      ) + (
            "cdef class gsl_ode_system(sage.gsl.ode.ode_system):\n"
	        if version < 2 else
	    "cdef class gsl_ode_system(sage.calculus.ode.ode_system):\n"
	      ) + ''.join(
            "    cdef Wrapper_rdf _flow%d\n"%i for i in il
              ) +
            "    def __init__( self, "
              + ', '.join(
	    "Wrapper_rdf f%d"%i for i in il
	      ) +
	    " ):\n"
              + ''.join(
            "        self._flow%d = f%d\n"% (i,i) for i in il
              ) +
            "    cdef int c_f( self, double t, double *y, double *dydt ):\n"
            #"        cdef double *result = [t]\n"
              + ''.join(
            "        self._flow%d.call_c( y, &dydt[%d] )\n"
            #"        dydt[%d] = result[0]\n"
                % (i,i) for i in il
              ) +
            "        return GSL_SUCCESS\n"
            # sometime include c_j() as well?
        )
        #print cython_code
        from sage.misc import cython
	try:
            module = sage.misc.cython.compile_and_load( cython_code )
	except RuntimeError as e:
	    if version > 0:
		return self.gsl_system( make_t_explicit=make_t_explicit, version=(version-1) )
	    else:
	        raise e
        T = ode_solver()
        #T.algorithm = 'bsimp'
        #print flow; sys.stdout.flush()
        print 'in gsl:\n', [ (v, flow[v]) for v in vars ]
        T.function = module.gsl_ode_system( *(
            fast_float( flow[v], *vars )
            for v in vars
        ) )
        return T
    def solve_gsl( self, initial_conditions, start_time=0, end_time=20, step=0.1, make_t_explicit=False ):
        ## don't cache the cython thing - it breaks load_session
        #try:
        #    self._gsl_system
        #except AttributeError:
        #    self._gsl_system = self.gsl_system()
        gsl_system = self.gsl_system( make_t_explicit=make_t_explicit )
        try:
            initial_conditions = [ initial_conditions(v) for v in self._vars ]
        except: pass
        if make_t_explicit:
            initial_conditions += [ start_time ]
        gsl_system.ode_solve(
            t_span=[start_time, end_time],
            num_points = (end_time - start_time) / step + 1,
            y_0 = initial_conditions
        )
        #print 'gsl solution', gsl_system.solution
        return Trajectory( self, ( [ t ] + ys for t,ys in gsl_system.solution ) )

# used by NumericalODESystem.solve()
# unlike the standard odeint(), this fails when the du/dt function
# raises an exception
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
        dynamics models to be fairly general, and able to produce specific
        flows for a varying number of "population indices".

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
        #print 'set_population_indices:', xi
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
