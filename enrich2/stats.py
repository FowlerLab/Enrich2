#  Copyright 2016 Alan F Rubin
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

from __future__ import absolute_import
import pandas as pd
import numpy as np
import scipy.stats as stats


def ztest_pair(df1, df2):
    """
    Given two data frames df1 and df2 each with columns 'score' and 'SE', 
    return a data frame containing both sets of scores and SE, now 
    labeled as 'score1', 'SE1', 'score2', and 'SE2', and the z statistic 
    and raw p-value for comparing each variant.
    
    """
    shared = df1.loc[:, ('score', 'SE')].merge(df2.loc[:, ('score', 'SE')], how='inner', left_index=True, right_index=True, suffixes=('1', '2'))
    shared['z'] = np.abs(shared['score1'] - shared['score2']) / np.sqrt(shared['SE1'] **2 + shared['SE2'] **2)
    shared['pvalue_raw'] = 2 * stats.norm.sf(shared['z'])
    return shared



def ztest_single(df1, score, SE):
    """
    Given a data frames df1 with columns 'score' and 'SE' and a score 
    and SE for a single variant, return a data frame containing the 
    z statistic and raw p-value for comparing each variant in the df
    to the given score and SE (usually wild type).
    
    """
    result = df1.loc[:, ('score', 'SE')]
    result['z'] = np.abs(result['score'] - score) / np.sqrt(result['SE'] **2 + SE **2)
    result['pvalue_raw'] = 2 * stats.norm.sf(result['z'])
    return result

