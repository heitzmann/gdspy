import re
#import sphinx.ext.autodoc

def myfactory(app, what, name, obj, options, lines):
    regex = re.compile('-+$')
    i = len(lines) - 3
    while i > 1:
        if regex.match(lines[i].strip()):
            lines[i] = ''
            lines[i-1] = '**' + lines[i-1] + '**'
            i -= 2
        i -= 1

def setup(app):
    app.connect('autodoc-process-docstring', myfactory)
