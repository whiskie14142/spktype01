A supporting module for jplephem to handle data type 1

Compute position and velocity of a celestial small body, from a NASA 
SPICE SPK ephemeris kernel file of data type 1 (Modified Difference Arrays).
http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

You can get SPK files for many solar system small bodies from HORIZONS 
system of NASA/JPL: http://ssd.jpl.nasa.gov/x/spk.html

Usage:
    from spktype01 import SPKType01
    kernel = SPKType01.open('path')
    position, velocity = kernel.compute_type01(center, target, jd)
    
    where:
        center - SPKID of central body (10 for Sol)
        target - SPKID of target body
        jd - time for computation (Julian date)

Author: Shushi Uetsuki (whiskie14142)
This module has been developed based on jplephem.spk and FORTRAN source 
of the SPICE Toolkit of NASA/JPL/NAIF.
jplephem : https://pypi.python.org/pypi/jplephem
SPICE Toolkit : http://naif.jpl.nasa.gov/naif/toolkit.html
