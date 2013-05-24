from pprint import PrettyPrinter

class SMSPrettyPrinter(PrettyPrinter):
    def format(self, object, context, maxlevels, level):
        if isinstance(object, float):
            return ('%.1f' % object), True, False
        else:
            return PrettyPrinter.format(self, object, context,
                                        maxlevels, level)


def wrap_onspace(text, width):
    """
        A word-wrap function that preserves existing line breaks
        and most spaces in the text. Expects that existing line
        breaks are posix newlines (\n).
        """
    return reduce(lambda line, word, width=width: '%s%s%s' %
                  (line,
                   ' \n'[(len(line[line.rfind('\n')+1:])
                          + len(word.split('\n',1)[0]
                                ) >= width)],
                   word),
                  text.split(' ')
                  )

import math
def wrap_always(text, width):
    """A simple word-wrap function that wraps text on exactly width characters.
       It doesn't split the text in words."""
    return '\n'.join([ text[width*i:width*(i+1)] \
                       for i in xrange(int(math.ceil(1.*len(text)/width))) ])


import re
def wrap(text, width):
    """Similar to wrap_onspace, but enforces the width constraint:
        words longer than width are split."""
    wordRegex = re.compile(r'\S{'+str(width)+r',}')
    return wrap_onspace(wordRegex.sub(lambda m: wrap_always(m.group(),width),text),width)
