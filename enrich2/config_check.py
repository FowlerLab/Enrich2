#  Copyright 2016-2017 Alan F Rubin
#
#  This file is part of Enrich2.
#
#  Enrich2 is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Enrich2 is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Enrich2.  If not, see <http://www.gnu.org/licenses/>.

"""
Functions for identifying the type of
:py:class:`~enrich2.storemanager.StoreManager` derived object associated with a
given configuration object (decoded from a JSON file as described `here
<https://docs.python.org/2/library/json.html>`_).

"""


def is_experiment(cfg):
    """
    Check if the given configuration object specifies an 
    :py:class:`~enrich2.experiment.Experiment`.

    Args:
        cfg (dict): decoded JSON object

    Returns:
        bool: True if `cfg` if specifies an 
        :py:class:`~enrich2.experiment.Experiment`, else False.

    """
    if 'conditions' in list(cfg.keys()):
        return True
    else:
        return False


def is_condition(cfg):
    """
    Check if the given configuration object specifies a 
    :py:class:`~enrich2.condition.Condition`.

    Args:
        cfg (dict): decoded JSON object

    Returns:
        bool: True if `cfg` if specifies a 
        :py:class:`~enrich2.condition.Condition`, else False.

    """
    if 'selections' in list(cfg.keys()):
        return True
    else:
        return False


def is_selection(cfg):
    """
    Check if the given configuration object specifies a 
    :py:class:`~enrich2.selection.Selection`.

    Args:
        cfg (dict): decoded JSON object

    Returns:
        bool: True if `cfg` if specifies a 
        :py:class:`~enrich2.selection.Selection`, else False.

    """
    if 'libraries' in list(cfg.keys()):
        return True
    else:
        return False


def is_seqlib(cfg):
    """
    Check if the given configuration object specifies a 
    :py:class:`~enrich2.seqlib.SeqLib` derived object.

    Args:
        cfg (dict): decoded JSON object

    Returns:
        bool: True if `cfg` if specifies a :py:class:`~enrich2.seqlib.SeqLib` 
        derived object, else False.

    """
    if 'fastq' in list(cfg.keys()) or 'identifiers' in list(cfg.keys()):
        return True
    else:
        return False


def seqlib_type(cfg):
    """
    Get the type of :py:class:`~enrich2.seqlib.SeqLib` derived object 
    specified by the configuration object.

    Args:
        cfg (dict): decoded JSON object

    Returns:
        str: The class name of the :py:class:`~seqlib.seqlib.SeqLib` derived 
        object specified by `cfg`.

    Raises:
        ValueError: If the class name cannot be determined.

    """
    if 'barcodes' in cfg:
        if 'map file' in cfg['barcodes']:
            if 'variants' in cfg and 'identifiers' in cfg:
                raise ValueError("Unable to determine SeqLib type.")
            elif 'variants' in cfg:
                return "BcvSeqLib"
            elif 'identifiers' in cfg:
                return "BcidSeqLib"
            else:
                raise ValueError("Unable to determine SeqLib type.")
        else:
            return "BarcodeSeqLib"
    elif 'overlap' in cfg and 'variants' in cfg:
        return "OverlapSeqLib"
    elif 'variants' in cfg:
        return "BasicSeqLib"
    elif 'identifiers' in cfg:
        return "IdOnlySeqLib"
    else:
        raise ValueError("Unable to determine SeqLib type for configuration "
                         "object.")


def element_type(cfg):
    """
    Get the type of :py:class:`~enrich2.storemanager.StoreManager` derived 
    object specified by the configuration object.

    Args:
        cfg (dict): decoded JSON object

    Returns:
        str: The class name of the 
        :py:class:`~enrich2.storemanager.StoreManager` derived object specified 
        by `cfg`.

    Raises:
        ValueError: If the class name cannot be determined.

    """
    if is_experiment(cfg):
        return "Experiment"
    elif is_condition(cfg):
        return "Condition"
    elif is_selection(cfg):
        return "Selection"
    elif is_seqlib(cfg):
        return seqlib_type(cfg)
    else:
        raise ValueError("Unable to determine type for configuration object.")
