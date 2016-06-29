from sage.all import *

class Bindings(dict):
    """Definitions for variables and functions.

    When you have an ODE system, for example, du/dt = a - bu^2,
    and you apply the bindings {a:1, b:2}, you get du/dt = 1 - 2u^2.
    The former is better for working with symbolically, but you need
    the latter to integrate and get a numeric trajectory.

    We also represent equilibria using bindings, like {hat(u): sqrt(2)/2}.
    In the future we may use them for the ODE's flow as well, e.g. {u: 1-2*u^2}.
    """
    def __init__(self, *args, **named_args):
        #print 'Bindings:__init__', ('args = ', args,
        #    ', named_args =', dict(**named_args))
        self._function_bindings = FunctionBindings()
        for a in args:
            try: # works if it's a FunctionBindings
                a.merge_into_function_bindings(self)
            except AttributeError:
                try:
                    self.update(a) # works if it's a dict or a Bindings
                    try: # and in case it's a bindings
                        self._function_bindings.update( a._function_bindings )
                    except AttributeError: # ok, it wasn't a Bindings
                        pass
                except ValueError:
                    print 'Unrecognized initializer for bindings:', a
        self.update(dict(**named_args))
        for k,v in self.items():
            kk = symbolic_expression(k)
            if (kk != k):
                del self[k]
            self[kk] = symbolic_expression(v)
	self.apply_to_self()
        #print ' =>', self
    def __repr__(self):
        frepr = repr(self._function_bindings)
        if frepr != '':
            frepr = ', ' + frepr
        return '{%s%s}' % (self.inner_repr(), frepr)
    def inner_repr(self):
        return ', '.join( '%s %s %s'%(k, '->', v) for k,v in self.items() ) 
    def latex_text(self):
        return '\\begin{align*}\n%s\n\\end{align*}' % self.latex_inner()
    def latex_inner(self):
	ltx = ' \\\\\n'.join( '  %s &\\to %s' % (latex(key), latex(val)) for key, val in self.items() )
        try:
            flx = self._function_bindings.latex_inner()
	    if ltx == '':
		return flx
	    if flx != '':
		return ltx + '\\\\\n' + flx
        except: pass
        return ltx
    def substitute(self, expr):
        """Apply the bindings to an expression"""
	# if it has a bind method, use that
	try:
	    return expr.bind( self )
	except AttributeError: pass
        # if it's a Bindings or dict, apply ourself to all the values
        try:
            return expr.__class__( { k:self.substitute(v) for k,v in expr.items() } )
        except (AttributeError, ValueError): pass
        # or if it's a vector or something, apply to all entries
        try:
            return expr.apply_map( lambda x : self.substitute( x ) )
        except AttributeError: pass
        # in case expr is a plain string
        expr = symbolic_expression(expr)
        expr = expr.substitute_expression(self)
        #print ' => ', expr
        return self._function_bindings.substitute(expr)
    def __call__(self, expr):
        """Apply the bindings to an expression"""
        return self.substitute(expr)
    def __deepcopy__(self, _dict):
        """I'm not getting the results I want from deepcopy(bindings) --
        it seems to garble the function bindings' values -- so here's my own"""
        other = Bindings()
        other.update( self )
        other._function_bindings = deepcopy( self._function_bindings, _dict )
        return other
    ## semantics of combining bindings.
    ## + should be a direct sum or something. I don't know whether it is - I
    ## think currently both + and bind are an odd relative of composition.
    ## merge is the same as +
    ## bind should be different.
    ## B1.bind( B2 ) is the same as B2( B1 ) and it means
    ## B1 but with B2 applied to its predicates (rhs's).
    def merge_in_place(self, *args, **named_args):
        other = Bindings(*args, **named_args)
	## apply other to all rhs of self, and import other's bindings
	## into self (except for things self already has bindings for)
        #for k, v in self.items():
        #    self[k] = other.substitute(v)
        try: # is it a FunctionBindings?
            other.merge_into_function_bindings(self)
        except AttributeError: # if not
            #self.update(other)
	    for k,j in other.iteritems():
		if k not in self:
		    self[k] = j
            other._function_bindings.merge_into_function_bindings(self)
        self.apply_to_self()
	#print 'result of bind_in_place: ', self
        return self
    def apply_to_self(self):
        """after merging bindings together, we have to apply the bindings to
        each other so that all substitutions get done in a single pass"""
	## TODO: fix this detecting-fixed-point shit
	strs = None
	while strs != str(self):
		strs = str(self)
		for k,v in self.items():
		    self[k] = self.substitute(v)
		self._function_bindings.apply_bindings(self)
	return self
    def bind_in_place(self, *args, **named_args):
	## @@ todo: this should be different from the merge or + operation?
        #print 'bind_in_place:', self, ',', args, ', ', named_args
        other = Bindings(*args, **named_args)
	## apply other to all rhs of self.
	## in certain special cases (to work with load_session() issue)
	## apply certain SR-variable replacements to lhs as well
	for k,j in self.iteritems():
	    self[k] = other(j)
	for k,j in other.iteritems():
	    if k.is_symbol() and j.is_symbol() and str(k) == str(j):
		if j in self:
		    self[k] = self[j]
		elif k in self:
		    self[j] = self[k]
	return self
    def merge(self, *args, **xargs):
        """Combine with another set of bindings.  We assume that self is the
        bindings that have already been applied, and the other bindings are
        being applied afterward.  Thus self's bindings take priority, if there's
        any potential conflict."""
	return deepcopy(self).merge_in_place( *args, **xargs )
    def bind(self, *args, **xargs):
        return deepcopy(self).bind_in_place( *args, **xargs )
    def __add__(self, other):
        return self.merge( other )

# see http://trac.sagemath.org/ticket/17553
# limit() and substitute_function() don't play well together.
# Once maxima returns a formal limit, it never gets re-evaluated
# even when the limit can be easily taken.  This re-evaluates them.
limop = limit( SR('f(x)'), x=0 ).operator()
def simplify_limits( expr ):
    #print 'before simplify_limits:', expr
    expr = expr.substitute_function( limop, maxima_calculus.sr_limit )
    #print 'after simplify_limits:', expr
    return expr

from sage.symbolic.function_factory import function
class FunctionBindings(Bindings):
    """I am annoyed that we have to substitute functions in a different way
    from other things.  I hope to find a better way, by formalizing the idea
    that parameters sometimes depend on other things."""
    def __init__(self, *args, **named_args):
        for a in list( args ) + [ named_args ]:
	    for k,v in a.items():
                # don't store as before-and-after functions due to
                # function-pickling woes
                # http://trac.sagemath.org/ticket/17558
                # store as (function name, argument list): return value
                # i.e. with no function objects stored
                try:
                    # if it's a FunctionBindings or other dict { (<name>,<args>): <expr> }
                    # it's tricky to exclude string-valued ks from this
                    def rterr(): raise TypeError
                    self.update( { (isinstance(k,tuple) and k or rterr()):v } )
                except (TypeError, ValueError):
                    try:
                        # if it's a dict { <fn or name>: <fn or expr> }, including
                        # the kind you would use with substitute_function(),
                        # transform it to { (<name>,<args>):<expr> }
			try:
			    self.update( { (str(k),v.arguments()):SR(v) } )
			except AttributeError:
                            self.update( { (str(k),SR(v).arguments()):SR(v) } )
		    except TypeError:
		        # this comes up if it's a Piecewise object
		        # not sure what else
		        self.update( { (str(k),()):v } )
                    except ValueError:
                        print 'Unrecognized initializer for bindings:', {k:v}
    def __repr__(self):
        # would like to use unicode arrow u'\u2192' but causes output codec error
        return ', '.join( '%s(%s) %s %s' % (key[0], ','.join( str(k) for k in key[1] ), '->', str(val)) for key,val in self.items() )
    def latex_inner(self):
        return '\\\\\n'.join( '  %s(%s) &\\to %s' % ( key[0], ','.join( latex(a) for a in key[1] ), latex( val ) ) for key, val in self.items() )
    def substitute(self, expr):
        """Apply the bindings to an expression"""
        try:
            return expr.apply_map( lambda x : self.substitute( x ) )
        except AttributeError: pass
        #print 'substitute ', expr,' + ', self
        for k,v in self.items():
            # evade confusion between function('f', nargs=1)
            # and function('f'), since they get interchanged by
            # maxima and pickle operations
            # http://trac.sagemath.org/ticket/17503
	    try: expr = expr.substitute_function( function(k[0]), v.function( *k[1] ) )
	    except AttributeError: pass
            try: expr = expr.substitute_function( function(k[0], nargs=len(k[1])), v.function( *k[1] ) )
	    except AttributeError: pass
            expr = expr.substitute_function(function(k[0]),v)
            expr = expr.substitute_function(function(k[0],nargs=len(k[1])),v)
            #print ' => ', expr
        return simplify_limits( expr )
    def __deepcopy__(self, _dict):
        return FunctionBindings( self )
    def merge(self, other={}, **args):
        return Bindings(self).bind(other, **args)
    def merge_into_function_bindings(self, other):
        """Merge this FunctionBindings into other bindings, in place.
	
        self is a FunctionBindings while other is a regular Bindings.
        this function provides duck typing because Bindings doesn't have it"""
        #print 'merge',self,'into',other
        other._function_bindings.update(self)
    def apply_bindings(self, bindings):
        for k, v in self.items():
            self[k] = bindings(v)

