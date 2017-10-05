"""Jump processes and their master equations, at least a subset of them"""
## these 2 lines seem to be needed for doctesting
## http://webcache.googleusercontent.com/search?q=cache:E3ULGhsSTZMJ:ask.sagemath.org/question/7998/trouble-using-import-and-doctest-together/+&cd=6&hl=en&ct=clnk&gl=us
import sys
sys.path.append(".")
##

from sage.all import *
import dynamicalsystems, hamiltonian, bindings, stochasticdynamics
from differenceequation import *

def mk_var( x, *args ):
    return SR.symbol( x+'_'+ '_'.join( str(a).replace('/','_') for a in args ),
        latex_name = x+'_{%s}' % (''.join(latex(a) for a in args)) )

class JumpProcess( stochasticdynamics.FiniteDimensionalStochasticDynamics ):
    def __init__(self, vars, transitions,
            time_variable=SR.symbol('t'),
            bindings=bindings.Bindings()):
        self._vars = vars
        self._transitions = transitions
        self._time_variable = time_variable
        self._bindings = bindings
    def bind( self, *args, **aargs ):
        return deepcopy( self ).bind_in_place( *args, **aargs )
    def bind_in_place( self, *args, **aargs ):
        b = bindings.Bindings( *args, **aargs )
        self._bindings.merge_in_place( b )
        self._transitions = [ (r,b(m)) for r,m in self._transitions ]
        return self
    def solve( self, initial_state, start_time=0, end_time=30, grain=1 ):
        self._grain = grain
        try: initial_state = [ initial_state(x) for x in self._vars ]
        except TypeError: pass
        return super( JumpProcess, self ).solve( initial_state, start_time=start_time, end_time=end_time )
    def update_time_and_state(self, t, x):
        return self.update_time_and_state_w_queue(t,x)
    def update_time_and_state_uniform(self, t, x):
        try:
            self._transitions_fast
        except AttributeError:
            self._transitions_fast = [
                (r,fast_callable(m, vars=self._vars, domain=float))
                for r,m in self._transitions
            ]
        ## compute rates of all transitions from the current state
        ## and the total rate
        transrates = []
        total_r = 0
        for r,mf in self._transitions_fast:
            er = mf(*x)
            if er > 0:
                total_r += er
                transrates.append( (r,mf,total_r) )
        ## pick one uniformly by rate
        r_pick = RR.random_element( 0, total_r )
        r,mf,tr = next( ( (r,mf,tr) for r,mf,tr in transrates if tr > r_pick ) )
        xnext = list(vector(x) + vector(r)*self._grain)
        ## exponential waiting time given total rate
        import numpy.random
        tnext = t + numpy.random.exponential( self._grain/total_r )
        print tnext; sys.stdout.flush()
        return ( tnext, xnext )
    def update_time_and_state_w_queue(self, t, x):
        ## construct a persistent store of events, rates, and conditional
        ## waiting times (i.e. waiting time if it's the next event to happen).
        ##  w = { i => (r,m,rate,waiting time) }
        import numpy.random
        def resample_hazard(i):
            b,m,_,__ = self._discrete_event_store[i]
            r = m(*x)
            h = numpy.random.exponential( 1 )
            w = (h/r if r > 0 else oo)
            self._discrete_event_store[i][2:4] = (r,w)
        try:
            self._discrete_event_store
        except AttributeError:
            self._discrete_event_store = {
                i : [r,fast_callable(m, vars=self._vars, domain=float),0,0]
                for i,(r,m) in enumerate(self._transitions)
            }
            for i in self._discrete_event_store.iterkeys():
                resample_hazard(i)
            ## record which events' rates are affected by which variables
            self._events_by_var = {
                v : set() for v in self._vars
            }
            for i,(r,m) in enumerate(self._transitions):
                for v in m.variables():
                    self._events_by_var[v].add(i)
        ## construct queue of events sorted by waiting time
        ##  q = [ i ]
        #discrete_event_queue = sorted(
        #    self._discrete_event_store.iterkeys(),
        #    key = lambda i: self._discrete_event_store[i][4]
        #)
        ### pop the first event, update by doing that transition
        #i = discrete_event_queue.pop(0)
        ## no point if queue gets re-sorted every time, just find the minimum
        i = min(
            self._discrete_event_store.iterkeys(),
            key=lambda i: self._discrete_event_store[i][3]
        )
        w = self._discrete_event_store[i][3]
        t = t + w
        to_update = set()
        for (ix,xi),inc,v in zip( enumerate(x), self._discrete_event_store[i][0], self._vars ):
            if inc != 0:
                x[ix] += inc
                to_update.update(self._events_by_var[v])
        resample_hazard(i)
        try: to_update.remove(i)
        except KeyError: pass
        ## update the affected rates and waits
        ##  w[i]
        def update_rate(i):
            b,m,r,w = self._discrete_event_store[i]
            rr = m(*x)
            w = (w*r/rr if r > 0 and rr > 0 else oo)
            self._discrete_event_store[i][2:4] = (rr,w)
        for i in self._discrete_event_store.iterkeys():
            self._discrete_event_store[i][3] -= w
            if i in to_update:
                update_rate(i)
        print t+w; sys.stdout.flush()
        return ( t + w, x )
    def deterministic_flow(self):
        flow = { v:0 for v in self._vars }
        for r,m in self._transitions:
            for i,ri in enumerate(r):
                try: flow[self._vars[i]] += ri*m
                except KeyError:
                    print 'KeyError on', i,ri
                    print self._vars[i]
                    pass
        return flow
    def deterministic_ode(self, time_variable=SR.symbol('t'), bindings=dynamicalsystems.Bindings()):
        return dynamicalsystems.ODESystem( 
            self.deterministic_flow(),
            self._vars, time_variable, bindings=bindings
        )
    def approximate_deterministic_difference_equation(self, step=1, time_variable=SR.symbol('t'), bindings=dynamicalsystems.Bindings()):
        flow = self.deterministic_flow()
        mapp = { v:(v + step*flow[v]) for v in flow.keys() }
        return DifferenceEquationSystem( mapp, self._vars, step=step, time_variable=time_variable, bindings=bindings )
    def stochastic_states( self, N ):
        from itertools import product
        if not all( sum(r) == 0 for r,_ in self._transitions ):
            # discretize the state space
            ss = [ vector(l) for l in product( *((i/N for i in range(N+1)) for s in self._vars) ) ] #if sum(l) <= 1 ]
        else:
            print 'Reducing dimension to', tuple(self._vars[:-1])
            ss = [ vector(l + (1 - sum(l),)) for l in product( *((i/N for i in range(N+1)) for s in self._vars[:-1]) ) if sum(l) <= 1 ]
        for s in ss: s.set_immutable()
        return ss
    def stochastic_state_binding_function( self ):
        return lambda s: dynamicalsystems.Bindings( { zv:zs for zv,zs in zip( self._vars, s ) } )
    def master_equations( self, N, p_name='p' ):
        ssts = self.stochastic_states( N )
        svs = [ mk_var( p_name, *sst ) for sst in ssts ]
        #print svs
        bind_rate_to_state = self.stochastic_state_binding_function()
        flow = { sv:0 for sv in svs }
        for sst in ssts:
            sv = mk_var( p_name, *sst )
            bv = bind_rate_to_state(sst)
            for r,m in self._transitions:
                mv = bv( m ) * sv
                if mv != 0:
                    tv = mk_var( p_name, *(sst+vector(r)/N) )
                    #print bv
                    #print sv, '=>', tv, ',', mv
                    flow[sv] -= mv
                    if tv not in flow: flow[tv] = 0
                    flow[tv] += mv
        return dynamicalsystems.ODESystem(
            flow,
            svs
        )
    def backward_equations( self, N, q_name='q' ):
        ssts = self.stochastic_states( N )
        bind_state = self.stochastic_state_binding_function()
        bflow = { mk_var( q_name, *s ) :
            sum( bind_state(s)(m) * ( mk_var( q_name, *(s+vector(r)/N) ) - mk_var(q_name,*s) )
            for r,m in self._transitions
            if bind_state(s)(m) != 0 )
            for s in ssts
        }
        qs = [ mk_var( q_name, *sst ) for sst in ssts ]
        return dynamicalsystems.ODESystem( bflow, qs )
    def generator_matrix( self, N, rate_ring=QQ ):
        ## matrix of instantaneous transition rates
        ## each row sums to 0
        ssts = self.stochastic_states(N)
        bind_state = self.stochastic_state_binding_function()
        dim = len(ssts)
        mtx = matrix( rate_ring, dim, dim, 0, sparse=True )
        s_idx = { s:i for i,s in enumerate(ssts) }
        for i,s in enumerate(ssts):
            bs = bind_state(s)
            for r,m in self._transitions:
                ms = bs(m)
                if ms != 0:
                    t = s + vector(r)/N
                    t.set_immutable()
                    mtx[i,s_idx[t]] += ms
                    mtx[i,i] -= ms
        return mtx
    def markov_transition_matrix( self, N, transition_ring=QQ ):
        ## matrix of transition probabilities
        ## each column sums to 1
        ## TODO: return a markov process, rt just a matrix?
        ssts = self.stochastic_states(N)
        bind_state = self.stochastic_state_binding_function()
        mtx = zero_matrix( transition_ring, len(ssts), sparse=True )
        s_idx = { s:i for i,s in enumerate(ssts) }
        for s in ssts:
            ps = [ ((s+vector(r)/N),bind_state(s)(m)) for r,m in self._transitions if bind_state(s)(m) != 0 ]
            if ps == []: ps = [ (s,1) ]
            pt = sum( rt for s1,rt in ps )
            for sr in ps: sr[0].set_immutable()
            for s1,rt in ps:
                mtx[s_idx[s1],s_idx[s]] += rt/pt
        return mtx
    def hamiltonian_system( self, p_vars=[], reduce=False, return_h=False ):
        if len(p_vars) == 0:
            p_vars = [ mk_var('p', v) for v in self._vars ]
        ## todo: this reduction repeats code with stochastic_states()?
        if reduce and all( sum(r) == 0 for r,_ in self._transitions ):
            print 'Reducing dimensions to', tuple(self._vars[:-1])
            reduce_bindings = dynamicalsystems.Bindings( { self._vars[-1]: 1 - sum(self._vars[:-1]) } )
            if len(p_vars) == len(self._vars):
                reduce_bindings.merge_in_place( { p_vars[-1] : 0 } )
                p_vars = p_vars[:-1]
        else:
            reduce=False
            reduce_bindings = dynamicalsystems.Bindings()
        vp = vector(p_vars)
        H = sum(
            reduce_bindings(m) * ( exp(
                vector(r[:len(vp)]).dot_product( vp )
            ) - 1 )
            for r,m in self._transitions
        )
        if return_h: return H
        x_vars = self._vars[:len(p_vars)]
        return hamiltonian.HamiltonianODE( H, x_vars, p_vars, bindings=reduce_bindings )
    def hamiltonian( self, p_vars=[], **aargs ):
        aaargs = { 'p_vars':p_vars, 'return_h':True }
        aaargs.update( aargs )
        return self.hamiltonian_system( **aaargs )
    ## note there is vestigial lagrangian code elsewhere
