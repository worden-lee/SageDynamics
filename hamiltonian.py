from sage.all import *
import dynamicalsystems

class HamiltonianODE(dynamicalsystems.ODESystem):
    def __init__(self, H, x_vars, p_vars, time_variable=SR.symbol('t'), bindings=dynamicalsystems.Bindings()):
        self._H = H
        self._configuration_vars = x_vars
        self._momentum_vars = p_vars
	import itertools
        super(HamiltonianODE,self).__init__(
            dict( { x:diff(H,p) for x,p in zip(x_vars,p_vars) },
                **{ p:-diff(H,x) for x,p in zip(x_vars,p_vars) }
            ),
            list( itertools.chain( x_vars, p_vars ) ),
            time_variable=time_variable,
            bindings=bindings
        )
        print 'time_var is', self._time_variable
    def bind_in_place(self, *largs, **dargs):
	b = dynamicalsystems.Bindings( *largs, **dargs )
	super(HamiltonianODE,self).bind_in_place( b )
	self._H = b(self._H)
    def equilibria( self, *ranges, **opts ):
	if opts.get( 'solve_numerically', False ):
	    return super(HamiltonianODE,self).equilibria( *ranges, **opts )
	constraints = opts.get( 'constraints', [] )
	## special solving procedure, taking into account the exponential
	## form in p.
	## note this may be applicable only to box hamiltonians.  Consider
	## making a specialization of this class for that.
	add_hats = self.add_hats()
	ui = dynamicalsystems.indexer('u')
	u_vars = [ ui[i] for i,p in enumerate( self._momentum_vars ) ]
	p_subs = { add_hats(p):log(u) for u,p in zip(u_vars, self._momentum_vars) }
	u_eqns = [ 0 == add_hats(rhs).subs(p_subs).canonicalize_radical() for rhs in self._flow.itervalues() ] + [ c.operator()( add_hats(self._bindings( c.lhs() ) ).subs(p_subs).canonicalize_radical(), add_hats( self._bindings( c.rhs() ) ).subs(p_subs).canonicalize_radical() ) for c in constraints ]
	import sys
	print u_eqns; sys.stdout.flush()
	u_sols = solve(
	    u_eqns,
	    [ add_hats(x) for x in self._configuration_vars + u_vars ],
	    solution_dict = True
	)
	u_subs = { u:exp(add_hats(p)) for u,p in zip(u_vars, self._momentum_vars) }
	return [ dict( {
		    add_hats(x): d[add_hats(x)].subs(u_subs).canonicalize_radical()
		    for x in self._configuration_vars
	        }, **{
		    add_hats(p): log(d[u]).subs(u_subs).canonicalize_radical()
		    for u,p in zip(u_vars, self._momentum_vars)
	        }
	    ) for d in u_sols
	] 
