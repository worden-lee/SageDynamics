# requires: $(SageDynamics)/dynamicalsystems.py $(SageDynamics)/latex_output.py
from sage.all import *
import os
import sys
sys.path.append( os.environ['SageDynamics'] )
from dynamicalsystems import *
from latex_output import *

from sage.symbolic.relation import solve
from sage.calculus.calculus import limit
from sage.misc.latex import latex

def limits(expr, lims):
    for x,y in lims.items():
        args = {str(x):y}
        #print 'limit', args, expr
        expr = limit(expr, **args)
        #print ' =>',expr
    return expr

class AdaptiveDynamicsException(DynamicsException):
    def __init__( self, message, latex_str=None ):
        DynamicsException.__init__( self, message, latex_str )

class AdaptiveDynamicsModel(ODESystem):
    """Given a population dynamics model in which dynamics of population
    sizes are driven by certain parameters, if we suppose that some of those
    parameters are subject to mutation, we can construct an adaptive dynamics
    model describing how selection will drive those parameters to evolve.

    This class infers the evolutionary dynamics from a population dynamics
    model given a distribution of mutations to parameters.

    Note it won't work as is for a lot of population models - among other
    things, it needs to be able to solve the population dynamics' equilibrium
    symbolically, and it assumes there is a unique equilibrium.  It also
    doesn't handle branching points, extinctions, etc."""
    def __init__(self, popdyn_model, phenotype_indexers,
        equilibrium=None,
        early_bindings=Bindings(), late_bindings=Bindings()):
        """Infer adaptive dynamics from population dynamics given an
        evolving phenotypic parameter that determines the variation in the
        relevant quantities in the population dynamics.

        popdyn_model: an object of class PopulationDynamicsModel
        phenotype_indexers: list of "indexers" for variables describing each
          population of creatures' character, which can change by a
          process of mutation and selection.
        equilibrium: Bindings mapping population dynamics state variables
          to their equilibrium values.  If not supplied, the model's
          stable nontrivial equilibrium is used.  If that's not unique,
          an exception is raised.  Note that trying to find stable
          equilibria can overwhelm the computer if parameters are unbound.
        early_bindings: Whereas things that are simple variables in
          the population dynamics, say r_i, an intrinsic growth rate, now
          need to be treated as functions of phenotypic quantities, say
          r_i(u_i), this argument is a dict mapping variable names to their new
          functional representations, for instance
          {'r_i':sage.symbolic.function_factory.function('r_i',u_i)}.
        late_bindings: In some cases you can just provide an early_bindings
          that expands the parameters to their ultimate values, but in other
          cases it's too complex for the equation solver to deal with those.
          So you can provide a symbolic expression like r_i(u_i) in
          early_bindings - that will be used during the adaptive dynamics
          calculations - and then expand that out to the value of r_i in
          late_bindings, and get the right answer without bogging the
          computer down."""
        self._debug_output = latex_output_base( write_to_string() )
        self._popdyn_model = popdyn_model
        self._phenotype_indexers = phenotype_indexers
        self._early_bindings = early_bindings
        self._late_bindings = late_bindings
        self.calculate_adaptive_dynamics()
        """Provide a mapping from e.g. SR.symbol('Xhat_i') to the
        population dynamics equilibrium value of X_i."""
        if equilibrium is None:
            equilibria = self._popdyn_model.stable_nontrivial_equilibria()
            # there should just be one nontrivial solution
            # otherwise, we would have to handle all this more carefully
            if len(equilibria) != 1:
                self._debug_output.write( "All equilibria of population dynamics:" )
                self._debug_output.write_block( self._popdyn_model.equilibria() )
                raise AdaptiveDynamicsException( 'Population Dynamics does not have a unique stable nontrivial equilibrium', self._debug_output._output._str )
            equilibrium = Bindings( equilibria[0] )
            #equilibrium += { hat(k):v for k,v in equilibrium.items() }
        self._equilibrium = equilibrium
        self._debug_output.write( "Equilibrium of the population dynamics:" )
        self._debug_output.write_block( self._equilibrium )
        self._flow = dict( (k, equilibrium(v)) for k,v in self._flow.items() )
        super(AdaptiveDynamicsModel, self).__init__( self._flow, list(self._flow.keys()),
          bindings=early_bindings + late_bindings + popdyn_model._bindings + self._equilibrium )
        # after the late bindings it might need going over
        # to resolve the limits
        #print 'and after:', self._vars[0], '=>', self._flow[self._vars[0]]
        self._debug_output.write( 'The adaptive dynamics comes out to' )
        self._debug_output.write_block( self )
        self._debug_output.close()
        print self._debug_output._output._str
    def calculate_adaptive_dynamics(self):
        # We get the invasion exponent for u_i by adding an extra population,
        # which will soon be asymptotically identical to u_i.
        #print 'ad popdyn:', self._popdyn_model
        self._extended_system = deepcopy(self._popdyn_model)
        self._extended_system.set_population_indices(self._popdyn_model._population_indices + ['i'])
        self._extended_system.bind_in_place( self._early_bindings )
        #self._debug_output.write_block( self._extended_system )
        #self._debug_output.write( "And with bindings: " )
        #self._debug_output.write_block(  self._popdyn_model._bindings )
        #self._extended_system = self._extended_system.bind(self._popdyn_model._bindings)
        self._debug_output.write( "Original pop-dyn system: " )
        self._debug_output.write_block( self._popdyn_model )
        self._debug_output.write( "Extended pop-dyn system: " )
        self._debug_output.write_block( self._extended_system )
        X_i = self._extended_system._population_indexer['i']
        dXi_dt = self._extended_system._flow[X_i]
        # Now: The invasion exponent for u_i is 1/X_i dX_i/dt
        I_i = dXi_dt / X_i
        self._debug_output.write( 'The invasion rate for population $i$ is:\n\\begin{align*}\n',
            '  \\mathscr I = \\frac{1}{%s}\\frac{d%s}{dt}' % (latex(X_i), latex(X_i)), ' &= ',
            latex( I_i ), '\n\\end{align*}\n' )
        #I_i = I_i.limit(X_i = 0)
        #self._debug_output.write( ' &= ', latex( I_i ), '\n\\end{align}\n' )
        # The invasion exponent we want is lim_{u_i->u, X_i->0}dI/du_i
        self._debug_output.write( 'And the invasion exponent is\n\\begin{align*}\n', 
            '\\\\\n'.join( '  \\frac{\\partial\\mathscr I}{\\partial %s} = \\lim_{%s,%s\\to0}%s' %
                (latex(u[j]), ','.join('%s\\to %s' % (latex(uu['i']),latex(uu[j])) for uu in self._phenotype_indexers),
                 latex(X_i), latex(diff(I_i, u['i']))) for j in self._popdyn_model._population_indices
                for u in self._phenotype_indexers ),
            '.\n\\end{align*}\n')
        dI_du = [diff(I_i,u['i']) for u in self._phenotype_indexers]
        #dI_du = [ equilibrium(dI_dui) for dI_dui in dI_du ]
        dI_du = [ (u[j], Xj, limits(dI_dui,
              dict((uu['i'],uu[j]) for uu in self._phenotype_indexers)))
            for j, Xj in zip(self._popdyn_model._population_indices,
              self._popdyn_model.population_vars())
            for u, dI_dui in zip(self._phenotype_indexers, dI_du) ]
        #print 'as u_i->u_*:\n', join( (" dI/d%s: %s" % (u_j, dI_duj)
        #    for u_j, X_j, dI_duj in dI_du), '\n')
        dI_du = [ (u_j, X_j, dI_duj.limit(X_i = 0)) for u_j, X_j, dI_duj in dI_du ]
        self._debug_output.write( 'Which comes out to\n\\begin{align*}\n', 
            '\\\\\n'.join( '  \\frac{\\partial\\mathscr I}{\\partial %s} = %s' %
                (latex(u_j), latex(dI_duj)) for u_j, X_j, dI_duj in dI_du ),
            '.\n\\end{align*}\n')
        # The adaptive dynamics system is du/dt = \gamma \hat{X}_i dI/du.
        # This is a general adaptive dynamics for the one-species,
        # one-resource Mac/Lev system - we'll supply a specific mapping
        # from u to b, m, c at solve time.
        self._vars = [ u_j for u_j, X_j, dI_duj in dI_du ]
        add_hats = self._popdyn_model.add_hats()
        self._S = dict( (u_j, add_hats(dI_duj)) for u_j, X_j, dI_duj in dI_du )
        self._debug_output.write( '\\[ ', 
            '\\mathbf S', latex( column_vector( [ u_j for u_j in self._S.keys() ] ) ), ' = ',
            latex( column_vector( [ dI_duj for dI_duj in self._S.values() ] ) ), ' \\]\n' )
        gamma = SR.var('gamma')
        self._flow = dict( (u_j, gamma*add_hats(X_j*dI_duj))
            for u_j, X_j, dI_duj in dI_du )

import numpy # do this now, not during compute_flow

class NumericalAdaptiveDynamicsModel( NumericalODESystem, AdaptiveDynamicsModel ):
    """Where AdaptiveDynamicsModel derives a closed expression for the
    adaptive dynamics using an explicit solution to the population dynamics'
    equilibrium condition, this class is appropriate when the population
    dynamics equilibrium needs to be solved numerically.  It sets up the
    apparatus for solving the equilibrium and plugging it into the AD
    equations, and then does it just-in-time at each step of the integration."""
    def __init__(self, popdyn_model, phenotype_indexers,
        equilibrium_function=None, time_variable=SR.symbol('t'),
        early_bindings=Bindings(), late_bindings=Bindings()):
        """Infer adaptive dynamics from population dynamics given an
        evolving phenotypic parameter that determines the variation in the
        relevant quantities in the population dynamics.

        popdyn_model: an object of class PopulationDynamicsModel
        phenotype_indexers: list of "indexers" for variables describing each
          population of creatures' character, which can change by a
          process of mutation and selection.
        equilibrium_function: function returning equilibrium bindings given
          bound model.
        early_bindings: Whereas things that are simple variables in
          the population dynamics, say r_i, an intrinsic growth rate, now
          need to be treated as functions of phenotypic quantities, say
          r_i(u_i), this argument is a dict mapping variable names to their new
          functional representations, for instance
          {'r_i':sage.symbolic.function_factory.function('r_i',u_i)}."""
        self._debug_output = latex_output_base( write_to_string() )
        self._popdyn_model = popdyn_model
        self._phenotype_indexers = phenotype_indexers
        self._early_bindings = early_bindings
        self._late_bindings = late_bindings
        self._equilibrium_function = equilibrium_function
        self.calculate_adaptive_dynamics()
        #print 'ad flow:', self._flow
        #print '+ bindings:', early_bindings + late_bindings + popdyn_model._bindings
        # unlike AdaptiveDynamicsModel, don't include the equilibrium here,
        # wait until compute_flow() to plug it in
        super(NumericalAdaptiveDynamicsModel, self).__init__(
          self._vars,
          time_variable = time_variable,
          flow = self._flow,
          bindings=early_bindings + late_bindings + popdyn_model._bindings )
        #print '--> flow:', self._flow
        self._debug_output.write( 'The adaptive dynamics comes out to' )
        self._debug_output.write_block( self )
        self._debug_output.close()
    def equilibrium_function( self, u_bindings ):
        if self._equilibrium_function is not None:
            return self._equilibrium_function( self, self._popdyn_model, u_bindings )
        answer = None
        for eq in self._popdyn_model.equilibria():
            eq = Bindings( { k:u_bindings(v) for k,v in eq.items() } )
            if any( eq(v) != 0 for v in self._popdyn_model.population_vars() ):
                if answer is not None:
                    raise AdaptiveDynamicsException( "Model has multiple nontrivial equilibria" )
                answer = eq
        if answer is None:
            raise AdaptiveDynamicsException( "Could not find interior equilibrium" )
	return Bindings( answer )
    def compute_flow( self, u, t ):
        #print 'compute_flow:', u, t
        u_bindings = Bindings( dict( zip( self._vars, u ) ) )
        eq = self.equilibrium_function( u_bindings )
        #print 'eq:', eq
        #print 'eq Xhat_0:', eq(hat('X_0'))
        if any( eq( hat(x) ) <= 0 for x in self._popdyn_model.population_vars() ):
             raise AdaptiveDynamicsException( 'Inviable population dynamics equilibrium in compute_flow' )
        #print 'flow 0:', self._flow[self._vars[0]]
        #print 'eq flow:', eq( self._flow[self._vars[0]] )
        #print 'u eq flow:', u_bindings( eq( self._flow[self._vars[0]] ) )
        #sys.stdout.flush()
        #print 'compute flow:', [ u_bindings( eq( self._flow[v] ) ) for v in self._vars ]
        #sys.stdout.flush()
        #print 'compute flow =', [ maxima( u_bindings( eq( self._flow[v] ) ) ) for v in self._vars ]
        #sys.stdout.flush()
        dudt = numpy.array( [ N( u_bindings( eq( self._flow[v] ) ) ) for v in self._vars ], float )
        #print 'compute flow', u, dudt
        #sys.stdout.flush()
        return dudt

# using this for the above equilibrium_function argument
# gives you adaptive dynamics using the unique interior equilibrium
# of the population dynamics (to be precise: equilibrium for which all
# population variables are nonzero.  It doesn't care about other state
# variables).  It will throw an exception if there isn't exactly one
# interior equilibrium.
def find_interior_equilibrium( ad, pd, u_bindings ):
    if not hasattr( find_interior_equilibrium, 'cache' ):
        find_interior_equilibrium.cache = {}
    if ad not in find_interior_equilibrium.cache:
        xs = [ hat(x) for x in pd.population_vars() ]
        #print 'create equilibrium function cache'
        #print 'ad early bindings', ad._early_bindings
        #print 'ad bindings', ad._bindings
        equilibria = [ { k:ad._bindings( v ) for k,v in eq.items() } for eq in pd.equilibria() ]
        #print 'equilibria', equilibria
        find_interior_equilibrium.cache[ad] = (equilibria, xs)
    equilibria, xs = find_interior_equilibrium.cache[ad]
    answer = None
    for eq in equilibria:
        if any( u_bindings(eq[x]) == 0 for x in xs ): continue
        if answer is not None:
            raise AdaptiveDynamicsException( 'Population dynamics has multiple interior equilibria' )
        answer = u_bindings( eq )
    if answer is not None:
        #print 'equilibrium:', answer
        return Bindings( answer )
    raise AdaptiveDynamicsException( 'Failed to find interior equilibrium' )

