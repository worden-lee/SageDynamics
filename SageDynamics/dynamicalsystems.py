"""Dynamical systems classes"""

from sage.all import *
from string import join
from sage.misc.latex import _latex_file_
latex.add_to_preamble('\\usepackage{amsmath}')

class ODE_trajectory(SageObject):
    """A trajectory of a system of ordinary differential equations.

    This class remembers the ODE that made it, and the symbolic variables
    involved in the ODE, so that the trajectories of complex symbolic
    expressions can be evaluated and plotted."""
    def __init__(self, system, timeseries):
        """system: the ODE_system object that provided the dynamics
        timeseries: a list of dictionaries, each providing a set of
            variable->number mappings, including the time variable.
        """
        self._system = system
        self._timeseries = timeseries
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
        if (xlabel == -1): xlabel = xexpr
        if (ylabel == -1): ylabel = yexpr
        xexpr = self._system._bindings(xexpr)
        yexpr = self._system._bindings(yexpr)
        #print 'after substitution: %s vs. %s' % (str(xexpr), str(yexpr))
        LP = list_plot(
          [ (xexpr.substitute(dic),yexpr.substitute(dic)) for dic in self._timeseries ],
          plotjoined = True, **args )
        try:
            xlabel.expand() # find out if it's an expression object
            xlabel = '$%s$' % latex(xlabel)
        except AttributeError: pass # if it's a string, leave it alone
        try:
            ylabel.expand()
            ylabel = '$%s$' % latex(ylabel)
        except AttributeError: pass
        LP.axes_labels( [xlabel,ylabel] )
        if (filename != ''):
            LP.save(filename)
        return LP
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
    the latter to integrate and get a numeric trajectory."""
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
        return '%s, %s' % (super(Bindings,self).__repr__(), self._function_bindings.__repr__())
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
        # in case it's a matrix or something
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
    def merge(self, other={}, **args):
        """Combine with another set of bindings.  We assume that self is the
        bindings that have already been applied, and the other bindings are
        being applied afterward.  Thus self's bindings take priority, if there's
        any potential conflict."""
        return deepcopy(self).bind(other, **args)
    def bind(self, other={}, **args):
        #print 'merge bindings:', self, ',', other, ', ', dict(**args)
        other = Bindings(other, **args)
        #for k, v in self.items():
        #    self[k] = other.substitute(v)
        try: # is it a FunctionBindings?
            other.merge_into_function_bindings(self)
        except AttributeError: # if not
            self.update(other)
            other._function_bindings.merge_into_function_bindings(self)
        self.apply_to_self()
        return self
    def apply_to_self(self):
        """after merging bindings together, we have to apply the bindings to
        each other so that all substitutions get done in a single pass"""
        for k,v in self.items():
            self[k] = self.substitute(v)
        self._function_bindings.apply_bindings(self)

from sage.symbolic.function_factory import function
class FunctionBindings(Bindings):
    """I am annoyed that we have to substitute functions in a different way
    from other things.  I hope to find a better way, by formalizing the idea
    that parameters sometimes depend on other things."""
    def __init__(self, *args, **named_args):
        for a in args:
            try:
                self.update(a) # works if it's a dict
            except ValueError:
                print 'Unrecognized initializer for bindings:', a
        self.update(dict(**named_args))
        for k,v in self.items():
            try:
                kk = function(k)
                self[kk] = v
                del self[k]
            except TypeError: pass # if k is already a function
    def __repr__(self):
        return super(Bindings,self).__repr__()
    def latex_inner(self):
        return '\\\\\n'.join( '  %s(%s) &\\to %s' % ( latex(key), ','.join( latex(a) for a in val.args() ), latex( val( *val.args() ) ) ) for key, val in self.items() )
    def substitute(self, expr):
        """Apply the bindings to an expression"""
        try:
            return expr.apply_map( lambda x : self.substitute( x ) )
        except AttributeError: pass
        #print 'substitute ', expr,' + ', self
        for k,v in self.items():
            expr = expr.substitute_function(k,v)
            #print ' => ', expr
        return expr
    def __deepcopy__(self, _dict):
        return FunctionBindings( self )
    def merge(self, other={}, **args):
        return self.wrap().bind(other, **args)
    def wrap(self):
        return Bindings(self)
    def merge_into_function_bindings(self, other):
        """To merge FunctionBindings and regular bindings together."""
        #print 'merge',self,'into',other
        other._function_bindings.update(self)
    def apply_bindings(self, bindings):
        for k, v in self.items():
            self[k] = bindings(v).function(*v.args())

class ODE_system(SageObject):
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
            and functions.  It is assumed that these substitutions have 
            been applied to the flow, and this is used to interpret symbolic
            expressions for plotting, etc.
        """
        #print 'ODE_system init', (flow, vars, time_variable, bindings)
        self._flow = flow
        self._vars = vars
        self._time_variable = time_variable
        self._bindings = bindings
    def __repr__(self):
        """Output the system as a system of differential equations"""
        return join( ('%s -> %s'%(v,self._flow[v])
            for v in self._vars), '\n')
    def __deepcopy__(self, _dict):
        """There seems to be a weird bug when you deepcopy a symbolic expression.
        So I'm trying to avoid that."""
        other = copy( self )
        other._flow = dict()
        other._flow.update( self._flow )
        other._bindings = deepcopy( self._bindings )
        return other
    def _latex_(self):
        """Output the system as a system of differential equations in LaTeX form"""
        return '\\begin{align*}\n' + '\\\\\n'.join( 
             r'\frac{d%s}{d%s} &= %s'%(latex(v),latex(self._time_variable),latex(self._flow[v]))
                for v in self._vars ) + '\n\\end{align*}'
    def write_latex(self, filename):
        """Output the system in LaTeX form to a file"""
        ltxout = open( filename, 'w' )
        ltxout.write( _latex_file_( self, title='' ) )
        ltxout.close()
    def time_variable(self):
        """Provide the independent variable, for instance for plotting"""
        return self._time_variable
    def bind_in_place(self, binding):
        """Apply bindings to my formulas without copying to a new object.

        See bind(), below."""
        # make binding into a Bindings if needed
        try:
            binding.substitute( 't' )
        except AttributeError:
             binding = Bindings( binding )
        self._flow = dict( (k, binding(v)) for k,v in self._flow.items() )
        self._bindings = self._bindings.merge(binding)
    def bind(self, bindings):
        """If you create a system with various symbolic parameters, like
        dy/dt = ay^2 + by + c, or something, you can't numerically 
        integrate it as is, you have to give values to parameters a, b, and c.
        This method gives you a system just like self, but with parameters
        bound to the values provided.

        bindings: a Bindings object recording values for some variables and/or
        functions.
        """
        bound = deepcopy( self )
        bound.bind_in_place( bindings )
        return bound
    def solve(self, initial_conditions, end_points=20, step=0.1):
        """Use a numerical solver to find a concrete trajectory of the system.

        initial_conditions: a list of values for the time variable and the
        state variables (in that order).

        This function returns the numerical solution, but it also stores it
        in the ODE_system object's state, so that subsequent calls to plot()
        refer to the most recent solution generated."""
        #print "desolve: %s, %s, %s, %s, %s" % (
        #  [self._flow[v] for v in self._vars],
        #  self._vars, initial_conditions, self._time_variable,
        #  end_points )
        soln = desolve_system_rk4(
          [self._flow[v] for v in self._vars],
          self._vars, initial_conditions, ivar=self._time_variable,
          end_points=end_points, step=step )
        # make that timeseries of lists into timeseries of dictionaries
        timeseries = [ dict( (k,v) for k,v in
          zip([self._time_variable]+self._vars,point) )
          for point in soln ]
        return ODE_trajectory(self, timeseries)
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
    def equilibria(self):
        equil_eqns = [ 0 == rhs for rhs in self._flow.values() ]
        return solve( equil_eqns, *self._vars, solution_dict=True )
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
    def stable_equilibria(self):
        return [ eq for eq in self.equilibria() if self.is_stable( self.jacobian_matrix( eq ) ) ]
    def is_stable(self, matrix):
        for e in matrix.eigenvalues():
            if real( e ) >= 0: return false
        return true
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
        jacobians_b = [ self.jacobian_matrix( eq ) for eq in equilibria_b ]
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

class PopulationDynamicsSystem(ODE_system):
    """Specialization of ODE_system for population dynamics systems"""
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
    # don't set _population_indices directly, call this, to keep the flow in sync
    def set_population_indices(self, xi):
        self._population_indices = xi
        self._flow = self.flow()
        self._vars = self._nonpop_vars + self.population_vars()
        #print "setting flow:", self._flow
        #self._flow = dict( (k,self._bindings(v)) for k,v in self._flow.items() )
        #print "  flow =>", self._flow
    # subclasses have to provide the flow
    def flow(self): pass
    def nontrivial_equilibria(self):
        return [ eq for eq in self.equilibria() if sum( eq[x] for x in self.population_vars() ) != 0 ]
    def stable_nontrivial_equilibria(self):
        return [ eq for eq in self.nontrivial_equilibria() if self.is_stable( self.jacobian_matrix( eq ) ) ]
    def plot_ZNGIs( self, xlims, ylims, vars=None, **args ):
        """Plot populations' zero net growth isoclines, using ODE_system's
        plot_isoclines()."""
        if vars is None:
            vars = self.population_vars()
        return self.plot_isoclines( xlims, ylims, [ (v,0) for v in vars ], **args )
