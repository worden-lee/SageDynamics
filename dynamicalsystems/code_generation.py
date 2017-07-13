import dynamicalsystems
import bindings

## emit code implementing ODESystem's dynamics in R

## TODO: use with(as.list(st,pa)) as in https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf
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

## todo: implementation of JumpProcess in R?
