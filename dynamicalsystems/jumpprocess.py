"""Jump processes and their master equations, at least a subset of them"""
from sage.all import *
from string import join
from sage.misc.latex import _latex_file_
import dynamicalsystems

def mk_var( x, *args ):
    return SR.symbol( x+'_'+ '_'.join( str(a).replace('/','_') for a in args ),
        latex_name = x+'_{%s}' % (''.join(latex(a) for a in args)) )

class JumpProcess:
	def __init__(self, vars, transitions):
		self._vars = vars
		self._transitions = transitions
		self._state_indices = { s:i for i,s in enumerate(self._vars) }
	def deterministic_ode(self, time_variable=SR.symbol('t'), bindings=dynamicalsystems.Bindings()):
		flow = { v:0 for v in self._vars }
		for r,m in self._transitions:
			for i,ri in enumerate(r):
				try: flow[self._vars[i]] += ri*m
				except KeyError:
					print 'KeyError on', i,ri
					print self._vars[i]
					pass
		return dynamicalsystems.ODESystem( 
			flow, self._vars, time_variable, bindings=bindings
		)
	def stochastic_states( self, N ):
		from itertools import product
	        if not all( sum(r) == 0 for r,_ in self._transitions ):
		    # discretize the state space
		    ss = [ vector(l) for l in product( *((i/N for i in range(N+1)) for s in self._vars) ) if sum(l) <= 1 ]
	        else:
		    print 'Reducing dimension to', tuple(self._vars[:-1])
		    ss = [ vector(l + (1 - sum(l),)) for l in product( *((i/N for i in range(N+1)) for s in self._vars[:-1]) ) if sum(l) <= 1 ]
	        for s in ss: s.set_immutable()
	        return ss
	def stochastic_state_binding_function( self ):
		return lambda s: dynamicalsystems.Bindings( { zv:zs for zv,zs in zip( self._vars, s ) } )
	def eqns( self, N, p_name='p' ):
		ssts = self.stochastic_states( N )
		bind_rate_to_state = self.stochastic_state_binding_function()
		flow = { mk_var( p_name, *sst ):0 for sst in ssts }
		for sst in ssts:
			sv = mk_var( p_name, *sst )
			bv = bind_rate_to_state(sst)
			for r,m in self._transitions:
				mv = bv( m ) * sv
				if mv != 0:
					tv = mk_var( p_name, *(sst+vector(r)/N) )
					print bv
					print sv, '=>', tv, ',', mv
					flow[sv] -= mv
					flow[tv] += mv
		return dynamicalsystems.ODESystem(
			flow,
			flow.keys()
		)
