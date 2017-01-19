from sage.all import *
from dynamicalsystems import *

from sage.symbolic.relation import solve
from sage.calculus.calculus import limit
from sage.misc.latex import latex

def limits(expr, lims):
    for x,y in lims.items():
        args = {str(x):y}
        #print 'limit', args, expr
        expr = limit(expr, **args)
        #print ' =>',expr
    # hacky workaround for maxima bug - if they don't get limited,
    # just substitute them.  This may produce wrong output.  We're already
    # getting wrong output, so this just makes it wrong less of the time.
    # Not ideal, but useful for now.
    #expr = expr.subs( lims )
    expr = simplify_limits( expr ) # nice fn from dynamicalsystems
    print 'simplify to', expr
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
        early_bindings=Bindings(), late_bindings=Bindings(),
        workaround_limits=False):
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
          One way to avoid that is by using model.symbolic_equilibria().
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
          computer down.
        workaround_limits: see http://trac.sagemath.org/ticket/17428
          The Maxima system sometimes loses track of limit operations
          and leaves temporary variables in final results. Enable this
          option to fix them by doing a straight substitute operation
          after the limit operation."""
        try:
                self._debug_output
        except AttributeError:
                self._debug_output = latex_output_base( write_to_string() )
        self._popdyn_model = popdyn_model
        self._phenotype_indexers = phenotype_indexers
        self._early_bindings = early_bindings
        self._late_bindings = late_bindings
        self._workaround_limits = workaround_limits
        self.calculate_adaptive_dynamics()
        """Provide a mapping from e.g. SR.symbol('Xhat_i') to the
        population dynamics equilibrium value of X_i."""
        if equilibrium is None:
            equilibria = self._popdyn_model.stable_nontrivial_equilibria()
            # there should just be one nontrivial solution
            # otherwise, we would have to handle all this more carefully
            if len(equilibria) != 1:
                self._debug_output.write( "All equilibria of population dynamics:" )
                self._debug_output.write( self._popdyn_model.equilibria() )
                raise AdaptiveDynamicsException( 'Population Dynamics does not have a unique stable nontrivial equilibrium', self._debug_output._output._str )
            equilibrium = Bindings( equilibria[0] )
            #equilibrium += { hat(k):v for k,v in equilibrium.items() }
        self._equilibrium = equilibrium
        self._debug_output.write( "Equilibrium of the population dynamics:" )
        self._debug_output.write( self._equilibrium )
        self._flow = { k : equilibrium(v) for k,v in self._flow.items() }
        super(AdaptiveDynamicsModel, self).__init__(
            self._flow,
            self._vars, #list(self._flow.keys()),
            bindings=early_bindings + late_bindings + popdyn_model._bindings + self._equilibrium
        )
        # after the late bindings it might need going over
        # to resolve the limits
        #print 'and after:', self._vars[0], '=>', self._flow[self._vars[0]]
        self._debug_output.write( 'The adaptive dynamics comes out to' )
        self._debug_output.write( self )
        self._debug_output.close()
        #print self._debug_output._output._str
    def calculate_adaptive_dynamics(self):
        # We get the invasion exponent for u_i by adding an extra population,
        # which will soon be asymptotically identical to u_i.
        self._vars = []
        self._flow = {}
        self._S = {}
        gamma = SR.var('gamma')
        #print 'ad popdyn:', self._popdyn_model
        for resident_index in self._popdyn_model._population_indices:
            X_r = self._popdyn_model._population_indexer[ resident_index ]
            # TODO: deepcopy not deep enough?
            extended_system = deepcopy(self._popdyn_model)
            #self._extended_system.set_population_indices(self._popdyn_model._population_indices + ['i'])
            #print 'just before mutate:', extended_system.__class__.__name__, extended_system
            mutant_index = self.mutate( extended_system, resident_index )
            #mutant_index = extended_system.mutate( resident_index )
            extended_system.bind_in_place( self._early_bindings )
            #self._debug_output.write_block( self._extended_system )
            #self._debug_output.write( "And with bindings: " )
            #self._debug_output.write_block(  self._popdyn_model._bindings )
            #self._extended_system = self._extended_system.bind(self._popdyn_model._bindings)
            self._debug_output.write( "Original pop-dyn system: " )
            self._debug_output.write( self._popdyn_model )
            self._debug_output.write( "Extended pop-dyn system: " )
            self._debug_output.write( extended_system )
            X_i = extended_system._population_indexer[ mutant_index ]
            dXi_dt = extended_system._flow[X_i]
            # Now: The invasion exponent for u_i is 1/X_i dX_i/dt
            I_i = dXi_dt / X_i
            c = ('The invasion rate for mutant population is:\n\\begin{align*}\n' +
                '  \\mathscr I = \\frac{1}{%s}\\frac{d%s}{dt}' % (latex(X_i), latex(X_i)) + ' &= ' +
                latex( I_i ) + '\n\\end{align*}\n')
            print c
            self._debug_output.write( c )
            #I_i = I_i.limit(X_i = 0)
            #self._debug_output.write( ' &= ', latex( I_i ), '\n\\end{align}\n' )
            # The invasion exponent we want is lim_{u_i->u, X_i->0}dI/du_i
            self._debug_output.write( 'And the invasion exponent is\n\\begin{align*}\n', 
                '\\\\\n'.join( '  \\frac{\\partial\\mathscr I}{\\partial %s} = \\lim_{%s,%s\\to0}%s' %
                    (latex(u[resident_index]), ','.join('%s\\to %s' % (latex(uu[mutant_index]),latex(uu[resident_index])) for uu in self._phenotype_indexers),
                     latex(X_i), latex(diff(I_i, u[mutant_index])))
                    for u in self._phenotype_indexers ),
                '.\n\\end{align*}\n')
            dI_du = [diff(I_i,u[mutant_index]) for u in self._phenotype_indexers]
            #dI_du = [ equilibrium(dI_dui) for dI_dui in dI_du ]
            limdict = { uu[mutant_index]:uu[resident_index] for uu in self._phenotype_indexers }
            #print 'limit of', dI_du[0]
            print 'limit as', limdict
            dI_du = [ ( u[resident_index], limits( dI_dui, limdict ) )
                for u, dI_dui in zip(self._phenotype_indexers, dI_du) ]
            if self._workaround_limits:
                ## first, if any limit operators left in there,
                ## replace the limit by a direct substitution
                dI_du = [ ( u, d.substitute_function(
                    sage.calculus.calculus._limit, 
                    lambda e, f, t: ( e.subs({f:t}) if f in limdict else e )
                ) ) for u,d in dI_du ]
                ## second, if any of the bound variables left without the
                ## limit operator, substitute
                dI_du = [ (u,d.subs(limdict)) for u,d in dI_du ]
            #print 'as u_i->u_*:\n', join( (" dI/d%s: %s" % (u_j, dI_duj)
            #    for u_j, X_j, dI_duj in dI_du), '\n')
            print 'after those limits:\n  ', '\n  '.join(str(i) for u_j, i in dI_du)
            #from sage.interfaces.maxima_lib import maxima_lib
            #print maxima_lib( dI_du[0][1] )
            dI_du = [ (u_j, limits(dI_duj, {X_i: 0})) for u_j, dI_duj in dI_du ]
            if self._workaround_limits:
                dI_du = [ (u,d.subs({X_i:0})) for u,d in dI_du ]
            self._debug_output.write( 'Which comes out to\n\\begin{align*}\n', 
                '\\\\\n'.join( '  \\frac{\\partial\\mathscr I}{\\partial %s} = %s' %
                    (latex(u_j), latex(dI_duj)) for u_j, dI_duj in dI_du ),
                '.\n\\end{align*}\n')
            #print 'after limits:\n  ', '\n  '.join(str(i) for u_j, i in dI_du)
            # The adaptive dynamics system is du/dt = \gamma \hat{X}_i dI/du.
            add_hats = self._popdyn_model.add_hats()
            ## todo: S should be one vector per population?
            ## careful handling here in case more than one population is
            ## adapting a given variable
            for u_j, dI_duj in dI_du:
                if u_j not in self._vars: self._vars += [ u_j ]
                self._S[u_j] = self._S.get(u_j,0) + add_hats(dI_duj) 
                self._flow[u_j] = self._flow.get(u_j,0) + gamma*add_hats(X_r*dI_duj)
        self._debug_output.write( '\\[ ', 
            '\\mathbf S', latex( column_vector( [ u_j for u_j in self._S.keys() ] ) ), ' = ',
            latex( column_vector( [ dI_duj for dI_duj in self._S.values() ] ) ), ' \\]\n' )
    def fake_population_index(self):
        # LotkaVolterraAdaptiveDynamics uses this, to do partial derivatives
        # separately from each other
        return self._popdyn_model.fake_population_index()
    def mutate( self, pops, i ):
        return pops.mutate(i)
    def solve(self, initial_conditions, **opts):
        continue_on_extinction = False #opts.pop( 'continue_on_extinction', False )
        ## initial_conditions might be a list, make sure it's a bindings
        try:
            initial_conditions = Bindings( { v:iv for v, iv in zip(self._vars, initial_conditions) } )
        except: pass
        # check that initial populations are positive
        if any( initial_conditions(self._bindings(x)) <= 0 for x in self._popdyn_model.equilibrium_vars() ):
            print ( 'Aborting solve: initial population sizes are non-positive:',
                str( [ x == N(initial_conditions(self._bindings(x))) for x in self._popdyn_model.equilibrium_vars() ] ) )
            raise AdaptiveDynamicsException(
                'Initial population sizes are non-positive: ' +
                str( [ x == N(initial_conditions(self._bindings(x))) for x in self._popdyn_model.equilibrium_vars() ] )
            )
        trajectory = super(AdaptiveDynamicsModel,self).solve(initial_conditions, **opts)
        # check that populations stay positive
        xs = [ self._bindings(x) for x in self._popdyn_model.equilibrium_vars() ]
        #print 'timeseries:', trajectory._timeseries
        print 'eq vars', self._popdyn_model.equilibrium_vars()
        print 'self bindings:', self._bindings
        print 'xs:', xs
        print 'first pt of timeseries:', trajectory._timeseries[0]
        print 'X at first pt of timeseries:', [ trajectory._timeseries[0](x) for x in xs ]
        for i in range(len(trajectory._timeseries)):
            exts = [ x for x in xs if N(trajectory._timeseries[i](x)) <= 0 ]
            if len( exts ) > 0:
                if not continue_on_extinction:
                    print 'Adaptive dynamics cut short by extinction at ' + str( trajectory._timeseries[i] )
                    #print '_timeseries[i]:', trajectory._timeseries[i]
                    print 'at which, with', { x:self._bindings(x) for x in self._popdyn_model.equilibrium_vars() }
                    print 'X values are', { x:trajectory._timeseries[i](self._bindings(x)) for x in self._popdyn_model.equilibrium_vars() }
                    tsi = Bindings( { k:N(v) for k,v in trajectory._timeseries[i].iteritems() } )
                    print 'with', tsi
                    print 'equivalently ', { x:tsi(self._bindings(x)) for x in self._popdyn_model.equilibrium_vars() }
                    print 'WTF', [ type(v) for v in trajectory._timeseries[i].values() ], '<-->', [ type(v) for v in tsi.values() ]
                    trajectory._timeseries = trajectory._timeseries[:i]
                    break
                else: # this needs some thought
                    pass #for ix in self._
        return trajectory
    def bind_in_place(self, *bindings, **args):
        b = Bindings( *bindings, **args )
        super(AdaptiveDynamicsModel, self).bind_in_place( b )
        self._popdyn_model.bind_in_place( b )
        for k,v in self._S.items():
                self._S[k] = self._bindings( v )
        return self

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
        ## KILL ME - this init doesn't call AdaptiveDynamicsModel's init
        self._workaround_limits = False
        self.calculate_adaptive_dynamics()
        print 'ad flow:', self._flow
        print '+ bindings:', early_bindings + late_bindings + popdyn_model._bindings
        # unlike AdaptiveDynamicsModel, don't include the equilibrium here,
        # wait until compute_flow() to plug it in
        super(NumericalAdaptiveDynamicsModel, self).__init__(
          self._vars,
          time_variable = time_variable,
          flow = self._flow,
          bindings=early_bindings + late_bindings + popdyn_model._bindings )
        print '--> flow:', self._flow
        self._debug_output.write( 'The adaptive dynamics comes out to' )
        self._debug_output.write( self )
        self._debug_output.close()
    def equilibrium_function( self, u_bindings ):
        if self._equilibrium_function is not None:
            return self._equilibrium_function( self, self._popdyn_model, u_bindings )
        return self.default_equilibrium_function( u_bindings )
    def default_equilibrium_function( self, u_bindings ):
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
    def compute_flow( self, u, t, with_exceptions=False ):
        #print 'ode:'
        #print self
        #print 'compute_flow:', u, t
        u_bindings = Bindings( dict( zip( self._vars, u ) ) )
        #print 'u_bindings:', u_bindings
        # note equilibrium_function can raise DynamicsExceptions if there's
        # a problem finding pop. dyn. equilibria
        eq = self.equilibrium_function( u_bindings )
        #print 'eq:', eq
        #print 'eq Xhat_0:', eq(hat('X_0'))
        if any( eq( hat(x) ) <= 0 for x in self._popdyn_model.population_vars() ):
            #print 'inviable equilibria', eq
            # to do: this is the place to detect one or more extinctions
            raise AdaptiveDynamicsException( 'Inviable population dynamics equilibrium ' + str(eq) + ' in compute_flow' )
        #print 'flow 0:', self._flow[self._vars[0]]
        #print 'eq flow:', eq( self._flow[self._vars[0]] )
        #print 'u eq flow:', [ u_bindings( eq( self._flow[v] ) ) for v in self._vars ]
        #sys.stdout.flush()
        import numpy
        try:
            dudt = numpy.array( [ N( u_bindings( eq( self._flow[v] ) ) ) for v in self._vars ], float )
        except TypeError:
            print 'trouble integrating ODE - can\'t evaluate', [ u_bindings( eq( self._flow[v] ) ) for v in self._vars ], 'numerically'
            raise
        #print 'compute flow', u, dudt
        sys.stdout.flush()
        if with_exceptions and self.detect_equilibrium( u, dudt, t ):
            # if there's an equilibrium, kick to the caller, who can
            # test for a branching point
            self._equilibrium_detected_at = t
            raise EquilibriumDetectedException()
        return dudt
    def solve(self, initial_conditions, start_time=0, end_time=20, step=0.1):
        self._equilibrium_detected_at = -oo
        try:
            traj = super(NumericalAdaptiveDynamicsModel,self).solve( initial_conditions, start_time, end_time, step )
        except TrajectoryInterruptedException, ex:
            print ex, ':', ex._reason
            if isinstance( ex._reason, EquilibriumDetectedException ):
                # in case we returned early due to equilibrium, see if
                # we're at an evolutionary branching point
                traj = ex._trajectory
                branched_state = self.test_for_branch( traj )
                if branched_state:
                    rest = self.solve( [ branched_state[v] for v in self._vars ], self._equilibrium_detected_at, end_time, step )
                    bt = traj + ODETrajectory( self, rest._timeseries )
                    #print 'composite trajectory after branching:', bt
                    return bt
                # if we're at equilibrium and didn't branch, skip ahead
                # to the end time
                end_point = dict( traj._timeseries[-1] )
                end_point[ self._time_variable ] = end_time
                traj.append( end_point )
                #print 'extended timeseries to', end_time
                #print 'end point is', end_point
            # TODO: if it returned early due to an extinction
            elif isinstance( ex._reason, UnboundedDynamicsException ):
                return ex._trajectory
            else:
                print 'unhandled TrajectoryInterruptedException:', ex._reason
                raise
        #print 'solve returns', traj._timeseries
        return traj
    def test_for_branch( self, traj ):
        # try branching various phenotypes, see if it goes
        import random
        indices = self._popdyn_model._population_indices
        random.shuffle( indices )
        indexers = self._phenotype_indexers
        random.shuffle( indexers )
        for i in indices:
            for u in indexers: 
                # construct extended model with perturbed copy of one type
                branched_popdyn = deepcopy(self._popdyn_model)
                sport = max( self._popdyn_model._population_indices ) + 1
                branched_popdyn.set_population_indices( self._popdyn_model._population_indices + [ sport ] )
                branched_ad = self.__class__( branched_popdyn, self._phenotype_indexers, self._equilibrium_function, self._time_variable, self._early_bindings, self._late_bindings )
                branched_state = dict( traj._timeseries[-1] )
                #print 'branched_state before branching:', branched_state
                for uu in self._phenotype_indexers:
                    branched_state[uu[sport]] = branched_state[uu[i]]
                branched_state[u[sport]] += 0.001
                branched_state[u[i]] -= 0.001
                #print 'branched_state after branching:', branched_state
                # test for whether they diverge
                try:
                    #print 'test for branch at state', branched_state
                    import numpy
                    flow = branched_ad.compute_flow( numpy.array( [ branched_state[v] for v in branched_ad._vars ], float ), self._equilibrium_detected_at )
                except AdaptiveDynamicsException: #nope
                    #print 'exception!'
                    break
                flow_dict = dict( zip( branched_ad._vars, flow ) )
                #print 'flow is', flow_dict
                divergence = (flow_dict[u[sport]] - flow_dict[u[i]])/(branched_state[u[sport]] - branched_state[u[i]])
                if divergence > 0:
                    # yes, we're branching
                    #print 'yes!'
                    self._popdyn_model = branched_ad._popdyn_model
                    self._vars = branched_ad._vars
                    self._S = branched_ad._S
                    self._flow = branched_ad._flow
                    self._add_hats = None
                    return branched_state
                #print 'no!'
        return false

# using this for the above equilibrium_function argument
# gives you adaptive dynamics using the unique interior equilibrium
# of the population dynamics (to be precise: equilibrium for which all
# population variables are nonzero.  It doesn't care about other state
# variables).  It will throw an exception if there isn't exactly one
# interior equilibrium.
def find_interior_equilibrium( ad, pd, u_bindings ):
    if not hasattr( find_interior_equilibrium, 'cache' ):
        find_interior_equilibrium.cache = {}
    key = str(ad._flow)
    if key not in find_interior_equilibrium.cache:
        xs = [ hat(x) for x in pd.population_vars() ]
        #print 'create equilibrium function cache'
        #print 'ad early bindings', ad._early_bindings
        #print 'ad bindings', ad._bindings
        equilibria = [ { k:ad._bindings( v ) for k,v in eq.items() } for eq in pd.equilibria() ]
        #print 'equilibria', equilibria
        find_interior_equilibrium.cache[key] = (equilibria, xs)
    equilibria, xs = find_interior_equilibrium.cache[key]
    answer = None
    for eq in equilibria:
        if any( u_bindings(eq[x]) <= 0 for x in xs ): continue
        if answer is not None:
            raise AdaptiveDynamicsException( 'Population dynamics has multiple interior equilibria' )
        answer = u_bindings( eq )
    if answer is None:
        raise AdaptiveDynamicsException( 'Failed to find interior equilibrium' )
    if any( u_bindings(eq[x]) >= 1e15 for x in xs ):
        raise UnboundedDynamicsException( 'Population Dynamics has become unbounded' )
    #print 'equilibrium:', answer
    return Bindings( answer )

