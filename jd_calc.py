from __future__ import division
from __future__ import print_function

import numpy as np

def jd_calc( year, month, day ):

    a = np.floor((month - 14) / 12.0)
    jd = np.floor((1461 * (year + 4800 + a)) / 4.0)
    jd += np.floor((367 * (month - 2 - 12 * a)) / 12.0)
    x = np.floor((year + 4900 + a) / 100.0)
    jd -= np.floor((3 * x) / 4.0)
    jd += day - 2432075.5  # was 32075; add 2400000.5

    jd -= 0.5  # 0 hours; above JD is for midday, switch to midnight.

    return jd # modified julian day
