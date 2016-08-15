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
        return self
    def extend_with_action_in_place(self, u=SR.symbol('u')):
        self._u_var = u
        self._vars += [u]
        self._flow[u] = self.lagrangian()
        return self
    def extend_with_action(self, u=SR.symbol('u')):
        return deepcopy(self).extend_with_action_in_place()
    def lagrangian(self):
        ## note this is the lagrangian in terms of x and p - it should
        ## really be in terms of x and \dot{x}, which is harder
        return sum( p*self._flow[x] for x,p in zip(self._configuration_vars,self._momentum_vars) ) - self._H
    def transformed_hamiltonian(self):
        ## transformation from x,p to w,s introduced by
        ## Hu, Phys.Rev.A, 36:12, 5782--5790
        try: self._H_hu
        except AttributeError:
            sx = lambda b,s: 's'+s
            ss = [ dynamicalsystems.xform_symbol( pi, sx, sx ) for pi in self._momentum_vars ]
            sub_dict = { pi:log(si) for pi,si in zip( self._momentum_vars, ss ) }
            wx = lambda b,s: 'w'+s
            ws = [ dynamicalsystems.xform_symbol(xi,wx,wx) for xi in self._configuration_vars ]
            if set( ws ) == set( [ws[0]] ):
                wx = lambda b,s: 'w_'+b+s
                ws = [ dynamicalsystems.xform_symbol(xi,wx,wx) for xi in self._configuration_vars ]
            sub_dict.update( {
                xi:-si * wi
                for xi,si,wi in zip(self._configuration_vars,ss,ws)
            } )
            self._hu_bindings = dynamicalsystems.Bindings( sub_dict )
            self._H_hu = self._hu_bindings( self._H ).canonicalize_radical()
        return self._H_hu
    def equilibria( self, *ranges, **opts ):
        if opts.get( 'solve_numerically', False ):
            return super(HamiltonianODE,self).equilibria( *ranges, **opts )
        constraints = opts.get( 'constraints', [] )
        ## special solving procedure, taking into account the exponential
        ## form in p.
        ## note this may be applicable only to box-model Hamiltonians.  Consider
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
