import dynamicalsystems

## emit code implementing ODESystem's dynamics in R

## TODO: use with(as.list(st,pa)) as in https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf
def R_ode_fn( self, name='odefn' ):
    ## this translation is naive: it only works if all the left and right
    ## hand sides of the Sage expressions involved, when printed in string
    ## form, yield valid R syntax.
    ## also may behave weirdly if any variable or parameter is called dxdt
    ## or something
    vars = sorted(self._vars, key=str)
    params = set.union( *[ set( f.variables() ) for f in self._flow.values() ] ) - set(vars)
    params = sorted(params, key=str)
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
    code += '    dxdt <- list(c(\n'
    code += ',\n'.join( '      ' + str(v) + ' = ' + str(self._flow[v]) for v in vars ) + '\n'
    code += '    ))\n'
    code += '    return(dxdt)\n'
    code += '  })\n'
    code += '}\n'
    return code

dynamicalsystems.ODESystem.R_ode_fn = R_ode_fn

## todo: implementation of JumpProcess in R?
