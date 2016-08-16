from sage.all import *

class Bindings:
    """Definitions for variables and functions.

    When you have an ODE system, for example, du/dt = a - bu^2,
    and you apply the bindings {a:1, b:2}, you get du/dt = 1 - 2u^2.
    The former is better for working with symbolically, but you need
    the latter to integrate and get a numeric trajectory.

    We also represent equilibria using bindings, like {hat(u): sqrt(2)/2}.
    In the future we may use them for the ODE's flow as well, e.g. {u: 1-2*u^2}.
    """
    def __init__(self, *args, **named_args):
        """
        Creates a Bindings object from some sequence of dictionaries,
        other Bindings, compatible iterators, and/or key=value pairs.

        INPUT:

        0 or more of the following:
        - dict whose keys are strings or Sage expressions,
          and values are Sage expressions.
          String keys are converted to Sage variable names.
        - Bindings.
        - key=value pairs, interpreted the same as key/value pairs
          in a dict argument.

        OUTPUT:

        - Bindings object whose action is to transform a Sage expression
          by making all the indicated replacements, in the order given.
        
        Default value is the identity Bindings, which replaces nothing.

        key=value pairs are interpreted as arguments given after any
        unnamed arguments.

        EXAMPLES:

        A simple Bindings that binds a variable to a scalar value.

        ::
            sage: var('y')
            y
            sage: b = Bindings( y=2 )
            sage: b( sin(x+y) )
            sin(x + 2)

        A trivial Bindings doesn't do anything.

        ::
            sage: b = Bindings()
            sage: b( sin(x+y) )
            sin(x + y)

        Transformations are composed in the order the arguments are given.

        ::
            sage: b = Bindings( { y:'z' }, Bindings( z='w' ), w=2 )
            sage: b( sin(x+y) )
            sin(x + 2)

        Bindings knows how to replace functions as well as variables
        (using substitute_function() rather than substitute()).
        ## TODO: detect function arguments automatically?

        ::
            sage: from sage.symbolic.function_factory import function
            sage: f = function('f')
            sage: b = Bindings( FunctionBindings( f=cos ) )
            sage: b
            { f(x) -> cos(x) }
            sage: b( f(x) )
            cos(x)

        Error checking

        ::
            sage: b = Bindings( None )
            Traceback (most recent call last):
            ...
            ValueError: Unrecognized value for Bindings: None

        Transformations supplied in sequence are applied in sequence.
        The sequence is all unnamed arguments, in order, followed by
        a single dictionary of all named arguments.

        ::
            sage: b = Bindings( {'a':'b'}, {'b':'c'} )
            sage: b('a')
            c

        The `+` operator and synonymous `merge` functions cause
        bindings to be applied in sequence.

        ::
            sage: b = Bindings( a='b' ) + Bindings( b='c' )
            sage: b('a')
            c

        Because the contents of a dictionary have no sequence,
        transformations supplied within a dictionary are applied in
        parallel, i.e. they don't interact with each other.

        ::
            sage: b = Bindings( { 'a':'b', 'b':'c' } )
            sage: b('a + 2*b')
            b + 2*c

        ::
            sage: b = Bindings( a='b', b='c' )
            sage: b('a + 2*b')
            b + 2*c

        The internal representation of the bindings has any sequential
        transformations pre-applied.

        ::
            sage: b = Bindings( a='b' ) + Bindings( b='c' )
            sage: b
            { b -> c, a -> c }

        """
        self._dict = {}
        self._function_bindings = FunctionBindings()
        self.merge_in_place( *args, **named_args )
    def merge_in_place(self, *args, **named_args):
        """Modify a Bindings by combining with another set of bindings.

        INPUT: some sequence of Bindings, dicts, and/or named arguments.

        OUTPUT: The original Bindings, modified such that its action
        is equivalent to the original Bindings followed by all the
        transformations given in the arguments.

        The arguments are used in sequence, in the same way as in
        `__init__`.

        EXAMPLES:

        The action of the merged Bindings is equivalent to applying the
        original Bindings and then the added transformations.

            ::
                sage: y = var('y')
                sage: b = Bindings( x = y )
                sage: b.merge_in_place( x=1, y=2 )
                { y -> 2, x -> 2 }
                sage: b
                { y -> 2, x -> 2 }

        This sequential application includes function bindings.
        Function bindings have to be applied after variable bindings,
        not in parallel, because of Sage's implementation, so they
        are not pre-applied to variable bindings.

            ::
                sage: p,q,r,s = var('p,q,r,s')
                sage: from sage.symbolic.function_factory import function
                sage: f = function('f')
                sage: b1 = Bindings( p=f(q) )
                sage: b2 = FunctionBindings( f = r.function(x) )
                sage: b3 = Bindings( r=s )
                sage: b = Bindings(); b.merge_in_place( b1, b2, b3 )
                { ... }
                sage: b
                { p -> f(q), r -> s, f(x) -> s }
                sage: b( p )
                s

        Test merging a Bindings:

            ::
                sage: b = Bindings( a='b' )
                sage: b.merge_in_place( Bindings( a=1, b=2 ) )
                { ... }
                sage: b
                { b -> 2, a -> 2 }

        Test merging a dictionary:

            ::
                sage: b = Bindings( a='b' )
                sage: b.merge_in_place( { 'a':1, 'b':2 } )
                { ... }
                sage: b
                { b -> 2, a -> 2 }

        Test merging a FunctionBindings:

            ::
                sage: from sage.symbolic.function_factory import function
                sage: f = function('f')
                sage: g = function('g')
                sage: b = Bindings( FunctionBindings( f = g(x).function(x) ) )
                sage: b.merge_in_place( FunctionBindings( f = SR(1).function(x), g = SR(2).function(x) ) )
                { ... }
                sage: b
                { f(x) -> 2, g(x) -> 2 }

        Function bindings are always evaluated after variable bindings.
        The semantics of ordered combinations of function and variable
        bindings that interact with each other are unclear and should
        probably be considered undefined until more work is done on that.
            ::
                sage: b = Bindings( { x:'f(x)' }, FunctionBindings( f = (x^2).function(x) ) )
                sage: b
                { x -> f(x), f(x) -> x^2 }
                sage: b(x)
                x^2
                sage: b = Bindings( FunctionBindings( f = (x^2).function(x) ), { x:'f(x)' } )
                sage: b
                { x -> x^2, f(x) -> x^4 }
                sage: b(x)
                x^2
                sage: b = Bindings( { x:'f(x)' }, FunctionBindings( f = SR('g(f(x))').function(x) ) )
                sage: b
                { x -> f(x), f(x) -> g(f(x)) }
                sage: b.merge_in_place( Bindings( { x:'g(x)' }, FunctionBindings( g = (x^2).function(x) ) ) )
                { ... }
                sage: b
                { x -> f(g(x)), f(x) -> f(x^2)^2, g(x) -> x^2 }
                sage: b(x)
                f(x^4)^2

        Recursive transformations are applied only as many times as given:

        ::
                sage: b = Bindings( x='f(x)' )
                sage: b.merge_in_place(b)
                { ... }
                sage: b
                { x -> f(f(x)) }

        SEEALSO:
            `+`
            `merge`
            `bind`
        """
        for a in list( args ) + [ named_args ]:
            try:
                ## case: a is a FunctionBindings
                ## 1. apply a to self's rhs's
                ## 2. add a to self's function bindings
                a.merge_into_function_bindings(self)
            except AttributeError:
                try:
                    try:
                        ## case: a is a Bindings
                        ## for each transformation in a:
                        ## apply it to self
                        ## add it to self

                        ## 1. apply variable bindings to self
                        for k,v in self._dict.iteritems():
                            self._dict[k] = v.substitute( a._dict )
                        for k,v in self._function_bindings.iteritems():
                            self._function_bindings[k] = v.substitute( a._dict )
                        ## 2. apply function bindings to self and add them
                        ## to self's function bindings
                        a._function_bindings.merge_into_function_bindings( self )
                        ## 3. add variable bindings to self
                        for k,v in a._dict.iteritems():
                            if k not in self._dict:
                                self._dict[k] = v
                    except AttributeError:
                        ## case: a is a dict
                        ## 0. promote contents of a to SR
                        ##    because function bindings come after variable
                        ##    bindings, pre-apply my fb to its vb
                        xa = { symbolic_expression(k):self._function_bindings(v) for k,v in a.iteritems() }
                        ## 1. apply a to self's rhs's
                        for k,v in self._dict.iteritems():
                            self._dict[k] = v.substitute(xa)
                        for k,v in self._function_bindings.iteritems():
                            self._function_bindings[k] = v.substitute(xa)
                        ## 2. add a's transformations to self
                        for k,v in xa.iteritems():
                            if k not in self._dict:
                                self._dict[k] = v
                except (AttributeError, ValueError, TypeError):
                    raise ValueError('Unrecognized value for Bindings: %s' % str(a))
        return self
    def __repr__(self):
        """
        String representation of Bindings.

        It is represented as a sequence of transformations, each rendered
        as a key-value pair separated by '->'.  Variable transformations are
        printed before function transformations.

        INPUT: None.

        OUTPUT: String representation.

        EXAMPLES:

            ::
                sage: b = Bindings( x=2 ) + FunctionBindings( f=tan )
                sage: b
                { x -> 2, f(x) -> tan(x) }
        """
        irepr = ', '.join( '%s %s %s'%(k, '->', v) for k,v in self._dict.items() ) 
        frepr = self._function_bindings._inner_repr()
        if irepr != '' and frepr != '':
            frepr = ', ' + frepr
        return '{ ' + irepr + frepr + ' }'
    def _latex_(self):
        """
        `\LaTeX` representation of Bindings.

        It is represented as a sequence of transformations, each rendered
        as a key-value pair separated by '`\to`'.  Variable transformations are
        printed before function transformations.

        INPUT: None.

        OUTPUT: String representation.

        EXAMPLES:

            ::
                sage: b = Bindings( x=2 ) + FunctionBindings( f=tan )
                sage: latex(b)
                \begin{align*}
                   x &\to 2 \\
                   f(x) &\to \tan\left(x\right)
                \end{align*}
        """
        return self.latex_text()
    def latex_text(self):
        return '\\begin{align*}\n%s\n\\end{align*}' % self.latex_inner()
    def latex_inner(self):
        ltx = ' \\\\\n'.join( '  %s &\\to %s' % (latex(key), latex(val)) for key, val in self._dict.items() )
        try:
            flx = self._function_bindings.latex_inner()
            if ltx == '':
                return flx
            if flx != '':
                return ltx + ' \\\\\n' + flx
        except: pass
        return ltx
    def substitute(self, expr):
        """Apply the bindings to an expression.

        INPUT:

        - a Sage expression (in the Symbolic Ring)
        - a string that evaluates to a Sage expression
        - elements of other rings might work, if they can be used with
          substitute() and substitute_function()

        OUTPUT: transformed Sage expression.

        EXAMPLES:

        Substitution in a Sage expression

        ::
            sage: alpha, beta, y = var( 'alpha beta y' )
            sage: f = function( 'f' )
            sage: ex = alpha*x + beta*f(y)
            sage: b = Bindings( alpha=2, beta=3 ) + FunctionBindings( f = tan )
            sage: b(ex)
            2*x + 3*tan(y)

        Substitution using a string

        ::
            sage: b = Bindings( alpha=3 )
            sage: b( 'alpha + x' )
            x + 3

        Use of Bindings to bind parameters in a dynamical system

        ::
            sage: import dynamicalsystems
            sage: alpha, beta, gamma, delta, y = var( 'alpha beta gamma delta y' )
            sage: ode = dynamicalsystems.ODESystem( { x: alpha*x + beta*y, y: gamma*x + delta*y }, [x,y] )
            sage: ode
            x -> alpha*x + beta*y
            y -> gamma*x + delta*y
            sage: b = Bindings( alpha = -1, beta = 1/2, gamma = 1/2, delta = -1 )
            sage: b(ode)
            x -> -x + 1/2*y
            y -> 1/2*x - y

        Application of one Bindings to another applies the bindings to
        the argument's right hand sides.  (See also bind and merge methods,
        for discussion of multiple ways of composing Bindings.)

        ::
            sage: numerical_bindings = Bindings( alpha=1, beta=2 )
            sage: functional_form = FunctionBindings( f = (alpha + tan(beta*x)).function(x) )
            sage: numerical_bindings( functional_form )
            { f(x) -> tan(2*x) + 1 }
        """
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
        try:
            # in case expr is a plain string
            expr = symbolic_expression(expr)
        except TypeError: # and if it can't become an SR expression, just pass it
            return expr
        expr = expr.substitute(self._dict)
        return self._function_bindings.substitute(expr)
    def __call__(self, expr):
        """Apply the bindings to an expression.

        SEEALSO:

            :meth:`substitute`
        """
        return self.substitute(expr)
    def __deepcopy__(self, _dc_dict):
        """I'm not getting the results I want from deepcopy(bindings) --
        it seems to garble the function bindings' values -- so here's my own.

        TODO: write a test clarifying what this is for
        """
        other = Bindings()
        other._dict.update( self._dict )
        other._function_bindings = deepcopy( self._function_bindings, _dc_dict )
        return other
    def merge(self, *args, **xargs):
        """Combine with another set of bindings.

        This operation is the same as the `+` operation applied to a
        sequence of Bindings.  It produces a Bindings whose action is
        equivalent to applying all the given Bindings in order.
        
        INPUT: any values that can be passed to the Bindings constructor.

        OUTPUT: a Bindings whose action is equivalent to the original
        Bindings followed by all the transformations given in the
        arguments.

        EXAMPLES:

            ::
                sage: y = var('y')
                sage: b = Bindings( x = y )
                sage: b.merge( x=1, y=2 )
                { y -> 2, x -> 2 }
                sage: b
                { x -> y }

        SEEALSO:
            `+`
            `merge_in_place`
            `bind`
        """
        return deepcopy(self).merge_in_place( *args, **xargs )
    def bind(self, *args, **xargs):
        """Apply transformations to right hand side values of a Bindings.

        INPUT: anything that can be provided when constructing a Bindings.

        OUTPUT: transformed Bindings

        EXAMPLES:

        Use of bind() to bind a name to a number:

            ::
                sage: y = var('y')
                sage: b = Bindings( x=y )
                sage: b
                { x -> y }
                sage: b.bind( y=2 )
                { x -> 2 }
                sage: b
                { x -> y }

        SEEALSO:
            `merge`
            `bind_in_place`
        """
        return deepcopy(self).bind_in_place( *args, **xargs )
    def bind_in_place(self, *args, **named_args):
        """Modify bindings by applying transformations to its right hand side
        values.

        INPUT: anything that can be provided when constructing a Bindings.

        OUTPUT: the same object, now modified by transforming its values.

        EXAMPLES:

        bind_in_place leaves the object changed:

            ::

                sage: y = var('y')
                sage: b = Bindings( x=y )
                sage: b
                { x -> y }
                sage: b.bind_in_place( y=2 )
                { x -> 2 }
                sage: b
                { x -> 2 }

        When replacing symbol names, bind_in_place (and therefore also `bind`
        and the () operator) also operates on left hand sides.
        This is important when fixing the latex name of a symbol,
        which is important due to some bugginess in save_session().

            ::
                sage: b = Bindings( xhat=RDF(0.3) )
                sage: latex(b)
                \begin{align*}
                  \mathit{xhat} &\to 0.3
                \end{align*}
                sage: xhat = SR.symbol('xhat', latex_name='\hat{x}')
                sage: b.bind_in_place( xhat=xhat )
                { ... }
                sage: latex(b)
                \begin{align*}
                  \hat{x} &\to 0.3
                \end{align*}

        SEEALSO:
            :trac:`17559`
        """
        other = Bindings(*args, **named_args)
        ## apply other to all rhs of self.
        ## in certain special cases (to work with load_session() issue)
        ## apply certain SR-variable replacements to lhs as well
        for k,j in self._dict.iteritems():
            self._dict[k] = other(j)
        for k,j in other._dict.iteritems():
            if k.is_symbol() and j.is_symbol() and str(k) == str(j):
                if j in self._dict:
                    self._dict[k] = self._dict[j]
                elif k in self._dict:
                    self._dict[j] = self._dict[k]
        return self
    def __add__(self, other):
        """Composition of Bindings.

        INPUT:
        - multiple Bindings objects
        - a Bindings object plus something that can be used to construct
          a Bindings

        OUTPUT: a Bindings object whose action is equivalent to applying
        the summand operations in order from left to right.

        EXAMPLE:
            ::
                sage: (Bindings(a=1) + Bindings(b=2))('a+b')
                3

        SEEALSO:
            `merge`
            `bind`
        """
        return self.merge( other )

# see http://trac.sagemath.org/ticket/17553
# limit() and substitute_function() don't play well together.
# Once maxima returns a formal limit, it never gets re-evaluated
# even when the limit can be easily taken.  This re-evaluates them.
limop = None
def simplify_limits( expr ):
    global limop
    if limop is None: 
        limop = limit( SR('anonfunc(x)'), x=0 ).operator()
    expr = expr.substitute_function( limop, maxima_calculus.sr_limit )
    return expr

from sage.symbolic.function_factory import function
class FunctionBindings(Bindings,dict):
    """I am annoyed that we have to substitute functions in a different way
    from other things.  I hope to find a better way, by formalizing the idea
    that parameters sometimes depend on other things.

    EXAMPLES:

    Test construction using an expression
        ::
            sage: b = FunctionBindings( f = x^2 )
            sage: b
            { f(x) -> x^2 }
            sage: b('f(y)')
            y^2

    Test construction using a bare function
        ::
            sage: b = FunctionBindings( f = cos )
            sage: b
            { f(x) -> cos(x) }
            sage: b('f(y)')
            cos(y)

    Test construction using piecewise, which has at some points in time
    caused trouble
        ::
            sage: b = FunctionBindings( f = piecewise( [[(-oo,0),0], [(0,oo),x]] ) )
            sage: b
            { f(x) -> piecewise(x|-->0 on (-oo, 0), x|-->x on (0, +oo); x) }
            sage: b('f(-1)')
            0

    Test copy construction
        ::
            sage: b = FunctionBindings( FunctionBindings( f = cos ) )
            sage: b
            { f(x) -> cos(x) }

    Note: FunctionBindings supports only symbolic functions, not lambdas
        ::
            sage: b = FunctionBindings( f = lambda x:0 )
            Traceback (most recent call last):
            ...
            AttributeError: 'function' object has no attribute 'number_of_arguments'

    Note also that FunctionBindings distinguishes its functions'
    argument signatures. The question of argument signatures is
    possibly more of an inconvenience than anything else, partly
    because Sage makes function names unique on a global level.
    FunctionBindings takes steps to get the function signature
    right when calling `substitute_function()`, to make sure the
    substitution happens correctly, even when functions have been
    garbled by saving and loading from pickles.
    There is unavoidable ambiguity and functions with the same name
    but different signatures should be avoided.
        ::
            sage: from sage.symbolic.function_factory import function
            sage: b = FunctionBindings( f = function('g') ) + FunctionBindings( f = (x^2).function(x) )
            sage: b
            { f() -> g(), f(x) -> x^2 }
            sage: f = function('f')
            sage: b(f(x))
            g()
            sage: b(f())
            g()
            sage: y = var('y')
            sage: b(f(x,y))
            g()
            sage: b = FunctionBindings( f = function('g',nargs=1) ) + FunctionBindings( f = (x^2).function(x) )
            sage: b
            { f(x) -> g(x) }
    """
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
                        ## this comes up when given a bare function as in 
                        ## FunctionBindings( f = cos )
                        #self.update( { (str(k),()):v } )
                        nn = v.number_of_arguments()
                        if nn == 1: args = (x,)
                        else:
                            args = tuple( SR.symbol() for i in range(nn) )
                        self.update( { (str(k),args):v(*args) } )
                    except ValueError:
                        raise TypeError( 'Unrecognized initializer for FunctionBindings:{}'.format({k:v}) )
    def __repr__(self):
        return '{ ' + self._inner_repr() + ' }'
    def _inner_repr(self):
        # would like to use unicode arrow u'\u2192' but causes output codec error
        return ', '.join( '%s(%s) -> %s' % (key[0], ','.join( str(k) for k in key[1] ), str(val)) for key,val in self.items() )
    def latex_inner(self):
        return '\\\\\n'.join( '  %s(%s) &\\to %s' % ( key[0], ','.join( latex(a) for a in key[1] ), latex( val ) ) for key, val in self.items() )
    def substitute(self, expr):
        """Apply the bindings to an expression"""
        try:
            return expr.apply_map( lambda x : self.substitute( x ) )
        except AttributeError: pass
        expr = symbolic_expression(expr)
        for k,v in self.iteritems():
            # evade confusion between function('f', nargs=1)
            # and function('f'), since they get interchanged by
            # maxima and pickle operations
            # http://trac.sagemath.org/ticket/17503
            try: expr = expr.substitute_function( function(k[0]), v.function( *k[1] ) )
            except AttributeError: pass
            try: expr = expr.substitute_function( function(k[0], nargs=len(k[1])), v.function( *k[1] ) )
            except AttributeError: pass
            #expr = expr.substitute_function(function(k[0]),v)
            #expr = expr.substitute_function(function(k[0],nargs=len(k[1])),v)
        return simplify_limits( expr )
    def __deepcopy__(self, _dict):
        return FunctionBindings( self )
    def merge(self, other={}, **args):
        return Bindings(self).bind(other, **args)
    def merge_into_function_bindings(self, subject):
        """Merge this FunctionBindings into other bindings, in place.
        
        self is a FunctionBindings while subject is a regular Bindings.
        this function provides duck typing because Bindings doesn't have it"""
        ## 1. apply self to subject's variable bindings
        ## Correction: Do not
        #subject._dict = { k:self(v) for k,v in subject._dict.iteritems() }
        ## 2. apply self to subject's function bindings
        ## 3. add self's bindings to subject's function bindings
        subject._function_bindings.merge_function_bindings_in_place(self)
    def merge_function_bindings_in_place( self, object ):
        """Merge other FunctionBindings into self, in place."""
        ## 1. apply object to self's bindings
        for tup,v in self.iteritems():
            self[tup] = object(v)
        ## 2. add object's bindings to self's
        for tup,v in object.iteritems():
            if tup not in self:
                self[tup] = v
    def __add__(self, other):
        return Bindings( self, other )

