from sage.all import *
from sage.misc.latex import *

# these are from sage.misc.latex
COMMON_HEADER = \
r'''\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{amsfonts}
\usepackage{graphicx}
\pagestyle{empty}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{latexml}
'''

LATEX_HEADER = (
r'''\documentclass{article}
''' + COMMON_HEADER +
r'''\oddsidemargin 0.0in
\evensidemargin 0.0in
\textwidth 6.45in
\topmargin 0.0in
\headheight 0.0in
\headsep 0.0in
\textheight 9.0in
''')

## two variants of latex(thing) allowing for text mode vs. math mode
## implementations
def latex_math( o ):
    if isinstance(o, basestring): return o
    try: return o.latex_math()
    except AttributeError: return latex( o )

def latex_text( o ):
    if isinstance(o, basestring): return o
    try: return o.latex_text()
    except AttributeError: return '$' + latex_math(o) + '$'

class latex_output_base(SageObject):
    '''latex_output: class that can collect latex output and do some formatting'''
    def __init__(self, _output):
        self._output = _output
    def write(self, *args):
        '''Output text directly.  Unlike file.write() (apparently), we support
        multiple arguments to write().'''
        ## note this used to use the latex() rather than latex_text(), which
        ## is a breaking change.  calling code that uses this in text mode will
        ## now get different and likely broken output.
        for a in args:
            self._output.write( latex_text( a ) )
        return self;
    def write_inline(self, *args):
        '''Output latex representation of each argument, inline in math mode'''
        for o in args:
            self.write( latex_math(o) )
        return self;
    def write_block(self, *args):
        '''Output latex representation of each argument, set apart in \\[ \\]'''
        ## note doesn't really work for multiple rows - use write_align
        ## or dgroup
        self.write( '\n\\[' )
        self.write( '\\\\\n'.join( latex_math(o) for o in args ) )
        self.write( '\\]\n' )
        return self
    def write_environment(self, envname, *stuff):
        return self.write(
            '\n\\begin{', envname, '}\n  ',
            '\\\\\n  '.join( latex_math(a) for a in stuff ),
            '\n\\end{', envname, '}\n'
        )
    def write_align( self, *stuff ):
        return self.write_environment( 'align*', *stuff )
    def write_equality(self, *args):
        return self.write( '\n\\[\n  ', ' = '.join( latex_math(a) for a in args ), '\n\\]\n' )
    def write_equation(self, *args):
        return self.write_equality( *args )
    def write_equality_aligned(self, *args):
        '''For convenience: write that one thing is equal to another (or more).
        This is because write_block( a == b ) often just writes "false"...'''
        return self.write_environment( 'align*',
            latex_math(args[0]) + ' &= ' + '\\\\\n    &= '.join( latex_math(a) for a in args[1:] ) );
    def close(self):
        self._output.close()
        return self;

def environment( envname, *stuff, **keywords ):
    outer_mode=keywords.pop( 'outer_mode', 'text')
    inner_mode=keywords.pop( 'inner_mode', 'math' )
    between_text = keywords.pop( 'between_text', '\\\\\n  ' )
    l_fn = (latex_math if inner_mode == 'math' else latex_text)
    return wrap_latex( 
        '\n\\begin{' + envname + '}\n  ' +
        between_text.join( l_fn(a) for a in stuff ) +
        '\n\\end{' + envname + '}\n',
        outer_mode
    )

def latex_block( *stuff ):
    return wrap_latex(
        '\n\\[\n  ' +
        '\n  '.join( latex_math(a) for a in stuff ) +
        '\n\\]\n',
        'text'
    )

def align_eqns( *stuff ):
    return environment( 'align*',
        wrap_latex( latex_math( stuff[0] ) + ' &= ' +
            '\\\\\n  &= '.join( latex_math( a ) for a in stuff[1:] ),
            'math'
        )
    )

def dgroup_eqns( *stuff ):
    ## write latex for a group of aligned equations using dgroup
    ## to do: refactor using dgroup() below, with some kind of
    ## equation class?
    return wrap_latex(
        '\\iflatexml\n' +
        ## use align* instead when using latexml
        latex_text( align_eqns( *stuff ) ) +
        '\\else\n' +
        latex_text(
            environment( 'dgroup*',
                environment( 'dmath*',
                    latex_math( stuff[0] ) + ' = ' + latex_math( stuff[1] ),
                    outer_mode='math'
                ),
                *[ environment( 'dmath*',
                    ' = ' + latex_math( s ),
                    outer_mode='math'
                ) for s in stuff[2:] ],
                between_text = ''
            )
        ) +
        '\\fi\n',
        'text'
    )

def dgroup( things, op=None ):
    ## latex for a list of relational expressions, aligned using dgroup
    ## general form for things is a nested list:
    ##  [ [ expr, (op, expr), (op, expr), ... ],
    ##    [ expr, ... ] ]
    ## allows full choice of ops.  simpler is
    ##  [ [ expr, expr, ... ], [ expr, expr, ... ], ... ], op='='
    ## (can we allow op to be a sage operator rather than a string?
    #  looks nontrivial with how latex is done inside ginac)
    if op is None:
        return wrap_latex(
            '\\iflatexml\n' +
            latex_text( environment( 'align*', things ) ) +
            '\\else\n' +
            latex_text( environment( 'dgroup*',
                *[ environment( 'dmath*', s, outer_mode='math') for s in things[2:] ]
            ) ) +
            '\\fi\n',
            'text'
        )
    else:
        return wrap_latex(
            '\\iflatexml\n' +
            latex_text( environment( 'align*',
                *( latex_math( l[0] ) + ' ' + op + ' ' +
                ''.join(
                        ('\\\\\n  &'+op+ ' ').join( latex_math( e ) for e in l[1:] ) )
                        for l in things )
                ) ) +
            '\\else' +
            latex_text( environment( 'dgroup*',
                *( environment( 'dmath*',
                    latex_math(l[0]) + ' ' + op + ' ' +
                    ('\\\\\n  '+op+' ').join( latex_math(e) for e in l[1:] ),
                    outer_mode='math' ) for l in things ),
                between_text='' ) ) +
            '\\fi\n',
            'text'
        )

class latex_output_file(latex_output_base):
    '''latex_output_file: class to write a latex file.
    Takes care of document header and footer; provides functions write_latex()
    (output latex representations of things inline) and write_block() (output
    latex of things in blocks) as well as regular write().'''
    def __init__(self, _file):
        '''Construct this object given a file object.
        It's better to use the function latex_output( filename ) rather than
        calling this constructor directly.'''
        super(latex_output_file,self).__init__(_file)
        self.write( LATEX_HEADER )
        self.write( latex_extra_preamble() )
        self.write( '\\begin{document}\n' )
    def close(self):
        '''Write closing latex commands and close file'''
        self.write( '\n\\end{document}\n' )
        super(latex_output_file,self).close()

class write_to_string:
    '''Class to use in place of a file to capture the output of latex_output_file'''
    def __init__(self):
        self._str = ''
    def write(self, arg):
        self._str += arg
    def close(self): pass
    def flush(self): pass

def latex_output( filename ):
    '''Create a latex_output_file object writing to the specified filename'''
    return latex_output_file( open( filename, 'w' ) )

class wrap_latex( SageObject ):
    '''A helper class for inserting literal latex code into a context that requires
    an object with a latex() method'''
    def __init__(self, text, mode='math'):
        self._str = text
        self._mode = mode
    def _latex_(self):
        return self._str
    def latex_math(self):
        if self._mode == 'math': return self._str
        else: raise TypeError, 'Object does not have a math-mode LaTeX representation'
    def latex_text(self):
        if self._mode == 'math': return '$'+self._str+'$'
        else: return self._str

def xform_symbol(v, xform_str, xform_latex):
    """little utility function to do things like add a hat or similar thing
    to a Sage variable, e.g. make X_i into \hat{X}_i."""
    try:
        # is it a symbol?
        if not v.is_symbol():
            # no, it's a more complex expression
            raise ValueError( str(v) + ' is not a symbol' )
    except AttributeError:
        # it's a string
        v = SR.symbol(v)
    import re
    name = str(v)
    mn = re.match( '^([^^_]*)(.*?)$', name )
    if mn: name = xform_str( mn.group(1), mn.group(2) )
    #name = re.sub( '^([^_]*)(.*?)$', r'\g<1>'+flair+'\g<2>', str(v) )
    latex_name = latex(v)
    #latex_name = re.sub( r'^([^_]*)(.*?)$', r'\\'+flair+'{\g<1>}\g<2>', latex(v) )
    ml = re.match( r'^([^^_]*)(.*?)$', latex_name )
    latex_name = xform_latex( ml.group(1), ml.group(2) )
    return SR.symbol( name, latex_name=latex_name )

def add_flair(v, flair):
    def xform_str( base, subs ): return base + flair + subs
    def xform_latex( base, subs ): return '\\' + flair + '{' + base + '}' + subs
    return xform_symbol(v, xform_str, xform_latex)

def hat(v):
    return add_flair(v, 'hat')

def dot(v):
    return add_flair(v, 'dot')

def scriptedsymbol( base, superscripts=(), subscripts=() ):
    try:
        base.expand() # see if it's an expression
    except AttributeError: # if not
        base = SR.symbol( base ) # it is now
    name, latex_name = str(base), '{%s}'%latex(base)
    import sys
    #print base, 'sub', subscripts, 'super', superscripts; sys.stdout.flush()
    if len(superscripts) > 0:
        name += '^' + '^'.join(str(s) for s in superscripts)
        if len(superscripts) > 1:
            latex_name += '^{%s}' % ''.join('{%s}'%latex(s) for s in superscripts)
        else:
            latex_name += '^{%s}' % latex(superscripts[0])
    if len(subscripts) > 0:
        name += '_' + '_'.join(str(s) for s in subscripts)
        if len(subscripts) > 1:
            latex_name += '_{%s}' % ''.join('{%s}'%latex(s) for s in subscripts)
        else:
            latex_name += '_{%s}' % latex(subscripts[0])
    return SR.symbol( name, latex_name=latex_name )

def subscriptedsymbol( base, *subscripts ):
    return scriptedsymbol( base, subscripts=subscripts )

def superscriptedsymbol( base, *superscripts ):
    return scriptedsymbol( base, superscripts=superscripts )

# TODO: integrate with stuff above
def write_tex_inline( vname, lname=None, fname=None, bindings=None ):
    if lname is None: lname = '\\'+str(vname)
    if fname is None: fname = str(vname)
    ltx = latex_output_base( open( fname+'.value.tex-inline', 'w' ) )
    if bindings is None:
        vval = vname
    else:
        vval = bindings(vname)
    vstr = str(vval)
    import re
    if re.search('\.\d*0$',vstr):
        print 'change', vstr, ':',
        vstr = str(RDF(vval))
    if re.search('\.\d{4,}',vstr):
        print 'reduce', vstr, ':',
        vstr = str(N(vval, digits=3))
    if re.search('\.0$', vstr):
        print 'truncate', vstr, ':'
        vstr = str(ZZ(vval))
    print vstr
    if lname != '':
        ltx.write( '$' + lname + ' = ' + vstr + '$' )
    else:
        ltx.write( '$' + vstr + '$' )
    ltx.close()

class column_vector(SageObject): # v.column() is missing?
    '''column_vector( x ) behaves just like vector( x ) but draws itself
    as a nice column in latex'''
    def __init__(self, *args):
        self.vector = vector( *args )
    def __call__(self, *args):
        return self.vector( *args )
    def __getitem__(self, i):
        return self.vector[i]
    def __iter__(self):
        return iter(self.vector)
    def __add__(self, other):
        return column_vector( self.vector.__add__( other.vector ) )
    def __getattr__(self, attr):
        at = getattr( self.vector, attr )
        if not callable( at ):
            return at
        def method_wrapper( *args, **named_args ):
            rv = at( *args, **named_args )
            if type(rv) == type( self.vector ):
                return column_vector(rv)
            return rv
        return method_wrapper
    def _latex_(self):
        return ( '\\left(\\begin{array}{c}\n%s\n\\end{array}\\right)' %
          '\\\\\n'.join( '  %s' % latex(e) for e in self.vector ) )
    def transpose(self):
        return self.vector

def visual_D( indx, f, f_args ):
    #print 'D', indx, f, f_args
    if indx == [0] and len(f_args) == 1:
        return SR.symbol( '%s_prime_%s' % (
                str(f),
                str( f_args[0] )
            ),
            latex_name = '{%s}\'(%s)' % (
                latex_math( f ),
                latex_math( f_args[0] )
            )
        )
    else:
        return SR.symbol( 'partial_%s_%s_%s' % (
                ''.join( str(i) for i in indx ),
                str(f),
                '_'.join( str(i) for i in f_args )
            ),
            latex_name='\partial_{%s}%s(%s)' % (
                ''.join( latex_math(i) for i in indx ),
                latex_math(f),
                ', '.join( latex_math(e) for e in f_args )
            )
        )

# useful parent class: expression converter that doesn't
# do anything
from sage.symbolic.expression_conversions import SubstituteFunction
class IdentityConverter(SubstituteFunction):
    def __init__(self):
        pass
    def composition(self, ex, operator):
        # override the parent class's function replacing step
        return operator(*map(self, ex.operands()))

class latex_partials_representation(IdentityConverter):
    def derivative( self, ex, operator ):
        #print ex, ' - operator is ', operator
        return visual_D( operator.parameter_set(), operator.function(), ex.operands() )

## Object to convert a symbolic expression to one that has the same
## latex representation, except that greek letters are sorted before
## roman ones in products.
class GreekFirstLatex(IdentityConverter):
    #from sage.symbolic.function_factory import function
    #gmul = function( 'times', print_latex=
    def arithmetic( self, ex, operator ):
	if operator == (2*SR.symbol('x')).operator():
	    ## too simple? sort factors so that things whose latex string
	    ## starts with '\\' are before the pure alphabetical ones.
	    ll = sorted( ex.operands(), key=lambda v: latex(v).replace('\\',' ') )
	    ## 1 ??
	    if len( ll ) > 1: ll = [ l for l in ll if l != 1 ]
	    ## -1 ?? this is weird, guess I don't understand signs in pynac
	    lm = [ l for l in ll if l == -1 ]
	    if lm:
		ll = [ l for l in ll if l != -1 ]
		ll[0] *= reduce( lambda x,y:x*y, lm, 1 )
	    ## don't know a way to enforce order of arguments to multiply
	    ## operator, so create a fake variable whose latex string is
	    ## the desired product.
	    ## thus the expression returned by this converter is suitable
	    ## only for printing in latex, not for doing math or anything
	    ## else with.
	    Msym = SR.symbol( 'M_{}'.format( ZZ.random_element(1e+10) ), latex_name=' '.join(latex(v) for v in ll) )
	    print latex(ex), ' ==> ', latex(Msym)
	    return Msym
	    ##
	    factors = set( ex.operands() )
	    greek_factors = set( [ v for v in factors if
		( latex(v) == '\\'+str(v) or latex(v) == '\\lambda' )
	    ] )
	    return wrap_latex( ' '.join(
		[ latex(v) for v in ( sorted( greek_factors ) + sorted( factors - greek_factors ) ) ]
	    ), 'math' )
	else:
	    return operator( *map(self, ex.operands()) )

GFL_memo = None
## return an expression whose latex representation has the greek letters
## sorted before the roman letters.
## this expression can not be used for calculations, only for latex.
def greek_first_latex_ex(ex):
    global GFL_memo
    if GFL_memo is None: GFL_memo = GreekFirstLatex()
    try: return GFL_memo(ex)
    except AttributeError: # if passed a non-expression
        return ex

## return latex representation of an expression, with the greek letters
## sorted before the roman ones.
def greek_first_latex(ex):
    return latex( greek_first_latex_ex(ex) )

