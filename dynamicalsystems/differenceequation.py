"""Difference equation"""
from sage.all import *
from dynamicalsystems import *
from latex_output import *
from bindings import *

# =============================================================================
# Difference Equation class for Sage
# =============================================================================

class DifferenceEquationSystem(SageObject):
    """A system of difference equations.

    This object holds a symbolic representation of the difference equations,
    and can use it in numerical integration as well as making it available
    for symbolic manipulations.
    """
    def __init__(self, map, vars, time_variable=SR.var('t'), step=1.,
            bindings=Bindings()):
        """Construct from the basic data.

        map: a dictionary with one entry for each state variable.
            Keys are variables, and values are symbolic expressions.
        vars: a list of the state variables.  This is conceptually
            equivalent to flow.keys(), but it's important to give them
            an ordering.  These are variables, not strings.
        time_variable: typically t, but could be something else.
        step: the time increment per step
        bindings: optional substitutions of values for symbolic variables
            and functions.  This is used to interpret symbolic
            expressions for plotting, etc.
        """
        self._map = map
        self._vars = vars
        self._time_variable = time_variable
        self._step = step
        self._bindings = bindings
        self._map = { k:bindings(v) for k,v in self._map.items() }
    def __repr__(self):
        """Output the system as a system of difference equations"""
        return join( ('%s -> %s'%(v,self._map[v])
            for v in self._vars), '\n')
    def _latex_(self):
        """Output the system as a system of difference equations in LaTeX form"""
        return self.latex_text() # not correct in math mode!
    def latex_text(self):
        return ('\\iflatexml\n' +
            '\\begin{align*}\n' + '\\\\\n'.join(
            r'%s &\mapsto %s'%(latex(v),latex(self._time_variable),latex(self._map[v]))
                for v in self._vars ) + '\n\\end{align*}\n'
            '\\else\n' +
            '\\begin{dgroup*}\n\\begin{dmath*}\n' + '\\end{dmath*}\n\\begin{dmath*}\n'.join(
             r'%s \mapsto %s'%(latex(v),latex(self._time_variable),latex(self._map[v]))
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
    def map( self, at=None ):
        if at is None: at = self._vars
        try:
            at( self._vars[0] )
        except TypeError: # if it's not a Bindings, make it one
            at = Bindings( zip( self._vars, at ) )
        return vector( [ at(self._map[v]) for v in self._vars ] )
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
        self._map = { k:self._bindings(v) for k,v in self._map.items() }
        return self
    def solve(self, initial_conditions, start_time=0, end_time=20, step=1.):
        """Construct a concrete trajectory of the system.

        initial_conditions: list or Bindings of initial values for the state variables"""
        ## if initial_conditions is a Bindings, make a list
        try: initial_conditions = [ initial_conditions(x) for x in self._vars ]
        except TypeError: 
            ## or if it's a dict, make a list
            try: initial_conditions = [ initial_conditions[x] for x in self._vars ]
            except TypeError: pass
        from numpy import arange
        state = initial_conditions
        timeseries = []
        for t in arange( start_time, end_time+step, step ):
            if t > start_time:
                state = self.update_state( state, t )
            timeseries.append( Bindings( { self._time_variable:t }, { v:x for v,x in zip(self._vars,state) } ) )
        return Trajectory(self, timeseries)
    def update_state( self, state, t ):
        """Implement the map by transforming current state into next state

        state: numerical vector, entries in same order as self._vars, or Bindings
        t: time at next state"""
        try: state(self._vars[0])
        except TypeError:
            state = Bindings( dict( zip( self._vars, state ) ), { self._time_variable:t } )
        return [ state( self._map[v] ) for v in self._vars ]
    ## to do: put various methods of ODESystem into a shared parent class for
    ## use by this class
