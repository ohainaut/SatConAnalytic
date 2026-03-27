#!/usr/bin/env python3
'''SatConAnalytic - utility helpers shared across the package.
'''

import json, os, sys, logging
import numpy as np

#----------------------------------------------------------------------------
class _Dict(dict):
    '''convenience: access dict as dict.element.element'''
    __getattr__= dict.__getitem__
    __setattr__= dict.__setitem__
    __delattr__= dict.__delitem__

#----------------------------------------------------------------------------
def read_conanJson(myFile):
    '''read a json file; try first local, then fallback to the package Data dir.
    return a dict with the content'''

    try:
        # Try to open the file as provided (local file or absolute path)
        with open(myFile) as infile:
            return json.load(infile, object_hook=_Dict)
    except FileNotFoundError:
        # If file not found, try to locate it in the same directory as this script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        fallback_path = os.path.join(script_dir, "Data", myFile)
        with open(fallback_path) as infile:
            return json.load(infile, object_hook=_Dict)


#------------------------------------------------------------------------------
def init_logger(log):
    '''Configure the satcon logger with a file handler (DEBUG) and a
    console handler (INFO).  Safe to call multiple times.'''

    log.setLevel(logging.DEBUG)
    log.propagate = False
    if log.handlers:
        return

    log_format = logging.Formatter('[%(levelname)-8s %(name)s/%(funcName)s] %(message)s')

    file_handler = logging.FileHandler('satCon.log')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(log_format)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(log_format)

    log.addHandler(file_handler)
    log.addHandler(console_handler)


#------------------------------------------------------------------------------
class NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy arrays and scalar types."""
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        return super().default(obj)

