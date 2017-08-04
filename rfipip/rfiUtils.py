# Copyright [2017] [SKA SA]
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import division, print_function
import numpy as np
from scipy import ndimage
from ConfigParser import SafeConfigParser
import os
from skimage import filters


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


def parse_frank_ini_file(ini_file='', required_sections=None):
    """
    (Modified parse_ini_file)
    Parse an ini file into a dictionary. No checking done at all.
    :param ini_file
    :param required_sections
    :return: a dictionary containing the configuration
    """
    if (ini_file == '') and ('FRANKENINI' in os.environ.keys()):
        ini_file = os.environ['FRANKENINI']
    if required_sections is None:
        required_sections = []
    parser = SafeConfigParser()
    files = parser.read(ini_file)
    if len(files) == 0:
        raise IOError('Could not read the config file, %s' % ini_file)
    for check_section in required_sections:
        if not parser.has_section(check_section):
            raise ValueError(
                'The config file does not seem to have the required %s section?' % (
                    check_section,))
    config = {}
    for section in parser.sections():
        config[section] = {}
        for items in parser.items(section):
            config[section][items[0]] = items[1]
    return config


def rfi_threshold(data_array):
    """
    Calculate threshold based on Yen algorithm
    :param data_array: numpy array
    :return: 
    """
    return filters.threshold_yen(data_array)
