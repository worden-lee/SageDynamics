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
    def __init__(self, map, vars, time_variable=SR.var('t'), step=1,
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
        y |-> ay^2 + by + c, or something, you can't numerically
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
    def solve(self, initial_conditions, start_time=0, end_time=20, bindings=Bindings()):
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
        t = start_time
        while t <= end_time:
            if t > start_time:
                state = self.update_state( state, t, bindings=bindings )
            timeseries.append( Bindings( { self._time_variable:t }, { v:x for v,x in zip(self._vars,state) } ) )
            t += self._step
        return Trajectory(self, timeseries)
    def cython_map_function(self, parameters=None, extra_bindings=Bindings()):
        if parameters is None:
            parameters = list( set().union( m.variables() for m in self._map.values() ).difference( self._vars ) )
        self._parameters = parameters
        SR_to_cython = SRCythonConverter()
        cython_code = (
            "from cpython cimport array\n" +
            "import array\n"
            "def update_state( array.array state, double t, array.array parameters ):\n" +
            ''.join(
            "    %s = state[%d]\n" %(str(v),i) for i,v in enumerate(self._vars)
            ) +
            ''.join(
            "    %s = parameters[%d]\n" %(str(p),i) for i,p in enumerate(self._parameters)
            ) +
            ''.join(
            "    %s = %s\n" %(k,SR_to_cython(b)) for k,b in extra_bindings._dict.iteritems()
            ) +
            ''.join(
            "    state[%d] = %s\n" %(i,SR_to_cython(self._map[v].collect_common_factors()))
                for i,v in enumerate(self._vars)
            ) +
            "    return 0\n"
        )
        print cython_code
        return cython_code
    def update_state_compiled( self, state, t, bindings=Bindings() ):
        try: self._cython_module
        except AttributeError:
            from sage.misc import cython
            self._cython_module = sage.misc.cython.compile_and_load( self.cython_map_function() )
            #print "--compiled"
        import array
        try:
            ## if it's a Bindings, unroll it into an array for the cython
            state = array.array( 'd', [ state(v) for v in self._vars ] )
        except TypeError:
            ## if not, hope it's an array already
            pass
        parameters = array.array( 'd', [ bindings(p) for p in self._parameters ] )
        #print "--call"
        try:
            ## use duck typing to catch it if it's not an array yet
            self._cython_module.update_state( state, t, parameters )
        except TypeError:
            ## convert and redo
            state = array.array( 'd', state )
            #print "--call"
            self._cython_module.update_state( state, t, parameters )
        #print "--called"
        return state
    def update_state_fast( self, state, t ):
        """Implement the map by transforming current state into next state

        state: numerical vector, entries in same order as self._vars, or Bindings
        t: time at next state"""
        try: 
            ## if state is in Bindings or dict form, flatten into list form
            state = [ state(v) for v in self._vars ]
        except TypeError:
            pass
        try: self._fast_map
        except AttributeError:
            ## compile the map into fast_callable form for fast calling
            self._fast_map = [
                fast_callable( self._map[v], vars=self._vars )
                for v in self._vars
            ]
        return [ vm(*state) for vm in self._fast_map ]
    def update_state_naive( self, state, t ):
        try: state(vars[0])
        except TypeError:
            state = Bindings( dict( zip( self._vars, state ) ), { self._time_variable:t } )
        return [ state( self._map[v] ) for v in self._vars ]
    def update_state( self, state, t, bindings=Bindings() ):
        return self.update_state_compiled( state, t, bindings=bindings )
    ## to do: put various methods of ODESystem into a shared parent class for
    ## use by this class

class SRCythonConverter(sage.symbolic.expression_conversions.ExpressionTreeWalker):
    def __init__(self, ex=SR(0)):
        super(SRCythonConverter,self).__init__(ex)
    def pyobject(self, ex, obj):
        return str(float(obj))
    def symbol(self, ex):
        return str(ex)
    def arithmetic(self, ex, operator):
        from sage.symbolic.operators import arithmetic_operators
        opnm = arithmetic_operators[operator]
        ## special case: ((a^-1)*b) ==> b/a
        if opnm == '*':
            def is_denom(f):
                try:
                    if arithmetic_operators[f.operator()] == '^':
                        base, pow = f.operands()
                        if pow == -1:
                            return true
                except: ## it doesn't have all those components
                    pass
                return false
            def partition( tf, it ):
                ## utility: partition iterator according to boolean criterion
                y, n = [], []
                for z in it:
                    if tf(z): y.append(z)
                    else:     n.append(z)
                return y,n
            denom, num = partition( is_denom, ex.operands() )
            if len(num) > 0:
                st = '(%s)' %('*'.join(self(z) for z in num))
            else:
                st = '1'
            if len(denom) > 0:
                st += '/(%s)' %('*'.join(self(1/z) for z in denom))
            return st
        ## default case: put the operator name between the operands
        if opnm == '^': opnm = '**'
        return opnm.join( '(%s)'%self(z) for z in ex.operands() )
    def composition(self, ex, operator):
        return str(operator) + '(' + ','.join(self(z) for z in ex.operands()) + ')'
