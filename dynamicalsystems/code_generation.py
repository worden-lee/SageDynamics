import dynamicalsystems
import bindings

## emit code implementing ODESystem's dynamics in R

def R_ode_fn( self, name='odefn', auxiliary_variables=bindings.Bindings() ):
    ## this translation is naive: it only works if all the left and right
    ## hand sides of the Sage expressions involved, when printed in string
    ## form, yield valid R syntax.
    ## also may behave weirdly if any variable or parameter is called dxdt
    ## or something
    vars = sorted(self._vars, key=str)
    #params = set.union( *[ set( f.variables() ) for f in self._flow.values() ] ) - set(vars)
    #params = sorted(params, key=str)
    code = '#!/usr/bin/R\n'
    code += name + ' <- function( t, state, params ) {\n'
    #code += '  ## state\n'
    #for v in vars:
    #    code += '  ' + str(v) + " <- st[['" + str(v) + "']]\n"
    #code += '  ## parameters\n'
    #for p in params:
    #    code += '  ' + str(p) + " <- pa[['" + str(p) + "']]\n"
    #code += '  ## the derivative\n'
    ## with(as.list(st,pa)) as in https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf
    code += '  with(as.list(c(state,params)), {\n'
    ## if there are auxiliary values in the ODE that are actually
    ## longer expressions in the state and parameters, expand them
    ## before use. For example, if the ODE includes N, and N is
    ## actually shorthand for S+I.
    for k,v in auxiliary_variables._dict.iteritems():
        ## note if there are FunctionBindings in the Bindings they'll be ignored
        code += '    ' + str(k) + ' <- ' + str(v) + '\n'
    code += '    dxdt <- c(\n'
    code += ',\n'.join( '      ' + str(v) + ' = ' + str(self._flow[v]) for v in vars ) + '\n'
    code += '    )\n'
    code += '    dxdt <- dxdt[names(state)]\n'
    code += '    return(list(dxdt))\n'
    code += '  })\n'
    code += '}\n'
    return code

dynamicalsystems.ODESystem.R_ode_fn = R_ode_fn

def R_sample_code( self, odefn_name='odefn' ):
    vars = sorted( self._vars, key=str )
    params = set.union( *[ set( f.variables() ) for f in self._flow.values() ] ) - set(vars)
    params = sorted(params, key=str)
    code =  '#!/usr/bin/R\n'
    code += 'library(deSolve)\n'
    code += '\n'
    code += 'initial.conditions <- c(\n'
    code += '  ' + ',\n  '.join( str(v) + ' = 1' for v in vars ) + '\n'
    code += ')\n'
    code += '\n'
    code += 'parameters <- c(\n'
    code += '  ' + ',\n  '.join( str(p) + ' = 1' for p in params ) + '\n'
    code += ')\n'
    code += '\n'
    code += 'times <- seq( 2017, 2021, by=1/10 )\n'
    code += 'trajectory <- data.frame(lsoda(initial.conditions, times, ' + odefn_name + ', parms=parameters))\n'
    return code

dynamicalsystems.ODESystem.R_sample_code = R_sample_code

## emit code implementing stochastic event simulation in R
## from a JumpProcess object
## Gillespie method, code adapted from Drake and Rohani, "Simulating
## Stochastic Epidemics"
## http://courses.ecology.uga.edu/ecol8520-spring2017/wp-content/uploads/sites/8/2017/01/stochastic-simulation.pdf

import jumpprocess

def R_jump_process_fn( self, name='simulate', stepfnname='onestep', auxiliary_variables=bindings.Bindings() ):
    ## this translation is naive: it only works if all the left and right
    ## hand sides of the Sage expressions involved, when printed in string
    ## form, yield valid R syntax.
    ## also may behave weirdly if any variable or parameter is called dxdt
    ## or something
    vars = sorted(self._vars, key=str)
    #params = set.union( *[ set( f.variables() ) for f in self._flow.values() ] ) - set(vars)
    #params = sorted(params, key=str)
    code = '#!/usr/bin/R\n'
    code += '\n'
    code += '# function ' + stepfnname + '() is used by ' + name + '(), below\n'
    code += stepfnname + ' <- function( t, state, params ) {\n'
    #code += '  ## state\n'
    #for v in vars:
    #    code += '  ' + str(v) + " <- st[['" + str(v) + "']]\n"
    #code += '  ## parameters\n'
    #for p in params:
    #    code += '  ' + str(p) + " <- pa[['" + str(p) + "']]\n"
    #code += '  ## the derivative\n'
    ## with(as.list(st,pa)) as in https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf
    code += '  with(as.list(c(state,params)), {\n'
    ## if there are auxiliary values in the ODE that are actually
    ## longer expressions in the state and parameters, expand them
    ## before use. For example, if the ODE includes N, and N is
    ## actually shorthand for S+I.
    for k,v in auxiliary_variables._dict.iteritems():
        ## note if there are FunctionBindings in the Bindings they'll be ignored
        code += '    ' + str(k) + ' <- ' + str(v) + '\n'
    code += '    rates <- c(\n'
    code += ',\n'.join( '      ' + str(m) for r,m in self._transitions ) + '\n'
    code += '    )\n'
    code += '    jumps <- matrix( c(\n'
    code += ',\n'.join( '      ' + ', '.join(str(ri) for ri in r) for r,m in self._transitions ) + ' ),\n'
    code += '      ncol=' + str(len(self._vars)) + ', byrow=TRUE\n'
    code += '    )\n'
    code += '    rates_c <- cumsum(rates)\n'
    code += '    rate_total <- rates_c[' + str(len(self._transitions)) + ']\n'
    code += '    if ( rate_total <= 0 ) { stop( "No transitions can be made" ) }\n'
    code += '    U <- runif(1)\n'
    code += '    j <- min(which(rates_c >= U*rate_total))\n'
    code += '    state.new <- state[c(' + ','.join( "'"+str(v)+"'" for v in self._vars ) + ')] + jumps[j,]\n'
    code += '    delta.t <- rexp(n=1, rate=rate_total)\n'
    code += '    return(list(delta.t=delta.t,state=state.new))\n'
    code += '  })\n'
    code += '}\n'
    code += '\n'
    code += '# function ' + name + '()\n'
    code += '# simulate stochastic model in time\n'
    code += name + ' <- function( nstep, initcond, params, t0=0 ) {\n'
    code += '    # nstep: number of jumps of the process to simulate\n'
    code += '    # initcond: starting state.\n'
    code += '    #  should be formatted like c( ' + '=0, '.join( [ str(v) for v in self._vars[0:2] ] + ['...'] ) + ' )\n'
    code += '    # params: values of parameters.\n'
    code += '    #  format should be same as for initcond\n'
    code += '    # t0: starting time. Time will be advanced by a stochastic waiting time at every step.\n'
    code += '  t <- t0\n'
    code += '  x <- initcond\n'
    code += '  output <- array( dim=c(nstep+1,' + str(len(self._vars)+1) + ') )\n'
    code += '  colnames(output) <- c( ' + ', '.join( [ '\'t\'' ] + [ '\''+str(v)+'\'' for v in self._vars ] ) + ' )\n'
    code += '  output[1,] <- c(t,x)\n'
    code += '  tryCatch( {\n'
    code += '      for (k in 1:nstep) {\n'
    code += '        l <- ' + stepfnname + '(t, x, params)\n'
    code += '        t <- t + l$delta.t\n'
    code += '        x <- l$state\n'
    code += '        output[k+1,] <- c(t,x)\n'
    code += '      }\n'
    code += '    },\n'
    code += '    error = function(e) print(e)\n'
    code += '  )\n'
    code += '  output\n'
    code += '}\n'
    code += '\n'
    code += '# ' + name + '() returns a matrix with columns t, ' + ', '.join( [ str(v) for v in self._vars[0:2] ] + [ '...' ] ) + '\n'
    code += '# It can be made into a data frame by for example as_tibble()\n'
    return code

jumpprocess.JumpProcess.R_jump_process_fn = R_jump_process_fn
