"""
.. module:: smsPrettyPrinter
   :synopsis: missing

.. moduleauthor:: missing <email@example.com>

"""

import math
import re
from pprint import PrettyPrinter


class SmsPrettyPrinter(PrettyPrinter):
    """
    An instance of this class represents a printing facility."
    
    """
    def format(self, entity, context, maxlevels, level):
        """
        TODO: write docstring
        
        """
        if isinstance(entity, float):
            return ('%.1f' % entity), True, False
        else:
            return PrettyPrinter.format(self, entity, context, maxlevels,
                                        level)


def wrapOnspace(text, width):
    """
    Preserves existing line breaks and most spaces in the text. Expects that
    existing line breaks are posix newlines (\\n).
    
    """
    return reduce(lambda line, word, width=width: '%s%s%s' %
                  (line, ' \n'[(len(line[line.rfind('\n') + 1:]) + \
                                len(word.split('\n', 1)[0]) >= width)], word),
                  text.split(' '))


def wrapAlways(text, width):
    """
    Wraps text on exactly width characters. It doesn't split the text in words.
    
    """
    return '\n'.join([text[width * i:width * (i + 1)] \
                      for i in xrange(int(math.ceil(1.*len(text) / width)))])


def wrap(text, width):
    """
    Similar to wrapOnspace, but enforces the width constraint: words longer
    than width are split.
    
    """
    wordRegex = re.compile(r'\S{' + str(width) + r',}')
    return wrapOnspace(wordRegex.sub(lambda m: \
                                    wrapAlways(m.group(), width), text), width)
