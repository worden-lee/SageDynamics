"""Stochastic dynamics classes"""
from sage.all import *
from string import join
from sage.misc.latex import _latex_file_
import dynamicalsystems

class FiniteDimensionalStochasticDynamics(SageObject):
    """Generic class for discrete time stochastic dynamics simulations.
    Will primarily be used for finite dimensional state space.  If infinite
    dimensional turns out to have different requirements, class
    hierarchy might change.
    TODO: introduce DiscreteTimeDynamics parent or something."""
    def __init__( self, vars, time_variable=SR.symbol('t'), bindings=dynamicalsystems.Bindings() ):
        self._vars = vars
        self._time_variable = time_variable
        self._bindings = bindings
    def solve( self, initial_state, start_time=0, end_time=30, update_fn=None ):
        t,x = ( start_time, initial_state )
        soln = []
        try:
            while t <= end_time:
                soln.append( tuple( [t] + list(x) ) )
                if update_fn is not None:
                        t,x = update_fn( t, x )
                else:
                        t,x = self.update_time_and_state( t, x )
        except StopIteration:
            # this happens when we hit an absorbing state
            pass
        #print 'soln:', soln
        return dynamicalsystems.Trajectory( self, soln )
    def update_time_and_state( self, t, x ):
        raise NotImplementedError, "Child class of FiniteDimensionalStochasticDynamics needs to implement update_time_and_state()"

class DiscreteMarkovProcess( FiniteDimensionalStochasticDynamics ):
    def update_time_and_state( self, t, x ):
        return ( t + 1, self.update_state( x ) )
    def update_state( self, x ):
        raise NotImplementedError, "Child class of DiscreteMarkovProcess needs to implement update_state()"

class MarkovMatrix( DiscreteMarkovProcess ):
    def __init__( self, vars, matrix, time_variable=SR.symbol('t') ):
        super( MarkovMatrix, self ).__init__( vars, time_variable )
        self._matrix = matrix
    def update_state( self, x ):
        p = RR.random_element(0,1)
        cum = 0
        # TODO: make this work.  Need conversion to/from states, matrix
        # indices
