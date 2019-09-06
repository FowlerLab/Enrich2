#  Copyright 2016-2019 Alan F Rubin
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

import numpy as np


def rml_estimator(y, sigma2i, iterations=50):
    """Implementation of the robust maximum likelihood estimator.

        ::

            @book{demidenko2013mixed,
              title={Mixed models: theory and applications with R},
              author={Demidenko, Eugene},
              year={2013},
              publisher={John Wiley \& Sons}
            }

    """
    w = 1 / sigma2i
    sw = np.sum(w, axis=0)
    beta0 = np.sum(y * w, axis=0) / sw
    sigma2ML = np.sum((y - np.mean(y, axis=0)) ** 2 / (len(beta0) - 1), axis=0)
    eps = np.zeros(beta0.shape)
    for _ in xrange(iterations):
        w = 1 / (sigma2i + sigma2ML)
        sw = np.sum(w, axis=0)
        sw2 = np.sum(w ** 2, axis=0)
        betaML = np.sum(y * w, axis=0) / sw
        sigma2ML_new = sigma2ML * np.sum(((y - betaML) ** 2) * (w ** 2),
                                         axis=0) / (sw - (sw2 / sw))
        eps = np.abs(sigma2ML - sigma2ML_new)
        sigma2ML = sigma2ML_new
    var_betaML = 1 / np.sum(1 / (sigma2i + sigma2ML), axis=0)
    return betaML, var_betaML, eps
