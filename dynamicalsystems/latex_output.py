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
	    self.latex(args[0]) + ' &= ' + '\\\\\n    &= '.join( latex_math(a) for a in args[1:] ) );
    def close(self):
        self._output.close()
        return self;

def environment( envname, *stuff, **keywords ):
    outer_mode=keywords.pop( 'outer_mode', 'text')
    inner_mode=keywords.pop( 'inner_mode', 'math' )
    l_fn = (latex_math if inner_mode == 'math' else latex_text)
    return wrap_latex( 
        '\n\\begin{' + envname + '}\n  ' +
        '\\\\\n  '.join( l_fn(a) for a in stuff ) +
        '\n\\end{' + envname + '}\n',
        outer_mode
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
	        ) for s in stuff[2:] ]
	    )
        ) +
	'\\fi\n',
	'text'
    )

def dgroup( *stuff ):
    ## latex for a list of objects, aligned using dgroup
    return wrap_latex(
	'\\iflatexml\n' +
	latex_text( environment( 'align*', stuff ) ) +
	'\\else\n' +
	latex_text( environment( 'dgroup*',
	    *[ environment( 'dmath*', s, outer_mode='math') for s in stuff[2:] ]
	) ) +
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


