from __future__ import division, print_function
import numpy as np
from scipy import ndimage


def replaceChar(text):
    for ch in ['\\',
               '`',
               '*',
               '{',
               '}',
               '[',
               ']',
               '(',
               ')',
               '>',
               ' ',
               '+',
               '-',
               '.',
               '!',
               '/',
               '\'']:
        text = text.replace(ch, "")
    return text


def open_blob(data, open=None):
    """
    
    :param data: 
    :param open: numpy array with shape
    :return: 
    """
    if open:
        op_struck = open
    else:
        op_struck = np.ones((2, 2))
    return ndimage.binary_opening(data, structure=op_struck)


def close_blob(data, close=None):
    """
    
    :param data: 
    :param close: numpy array with shape
    :return: 
    """
    if close:
        cl_struck = close
    else:
        cl_struck = np.ones((1, 1))
    close_img = ndimage.binary_closing(data, structure=cl_struck)
    return np.invert(close_img)


def percentage(part, whole):
    """
    
    :param part: 
    :param whole: 
    :return: 
    """
    return 100 * float(part)/float(whole)
