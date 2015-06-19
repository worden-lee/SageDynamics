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

class latex_output_base(SageObject):
    '''latex_output: class that can collect latex output and do some formatting'''
    def __init__(self, _output):
        self._output = _output
    def latex(self, a):
        if isinstance(a, basestring):
            return a
        return latex(a)
    def write(self, *args):
        '''Output text directly.  Unlike file.write() (apparently), we support
        multiple arguments to write().'''
        for a in args:
            self._output.write( self.latex( a ) )
        return self;
    def write_inline(self, *args):
        '''Output latex representation of each argument, inline in math mode'''
        for o in args:
            self.write( '$%s$' % self.latex(o) )
        return self;
    def write_block(self, *args):
        '''Output latex representation of each argument, set apart in \\[ \\]'''
        self.write( '\n\\[' )
        self.write( '\\\\\n'.join( self.latex(o) for o in args ) )
        self.write( '\\]\n' )
        return self
    def write_environment(self, envname, *stuff):
        return self.write(
            '\n\\begin{', envname, '}\n  ',
            '\\\\\n  '.join( self.latex(a) for a in stuff ),
            '\n\\end{', envname, '}\n'
        )
    def write_align( self, *stuff ):
        return self.write_environment( 'align*', *stuff )
    def write_equality(self, *args):
        return self.write( '\n\\[\n  ', ' = '.join( self.latex(a) for a in args ), '\n\\]\n' )
    def write_equality_aligned(self, *args):
        '''For convenience: write that one thing is equal to another (or more).
        This is because write_block( a == b ) often just writes "false"...'''
        return self.write_environment( 'align*',
	    self.latex(args[0]) + ' &= ' + '\\\\\n    &= '.join( self.latex(a) for a in args[1:] ) );
    def close(self):
        self._output.close()
        return self;

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
    # maybe not needed any more?
    def __init__(self, text):
        self._str = text
    def _latex_(self):
        return self._str

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


