# spktype01
A supporting module for jplephem to handle data type 1 (Version 1.0.0)

This module computes positions and velocitys of a celestial small body, from a NASA SPICE SPK ephemeris kernel file of data type 1 (Modified Difference Arrays).
http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

You can get SPK files for many solar system small bodies from HORIZONS system of NASA/JPL.  See https://ssd.jpl.nasa.gov/?horizons

This module reads SPK files of data type 1, one of the types of binary SPK file.  

At the point of Oct. 2018, HORIZONS system provides files of type 21 for binary SPK files by default.  You can get type 1 binary SPK file for celestial small bodies through TELNET interface by answering back '1' for 'SPK file format'.

## Module required
* jplephem (version 2.6 or later)
* numpy

    Usage:
        from spktype01 import SPKType01
        kernel = SPKType01.open('path')
        position, velocity = kernel.compute_type01(center, target, jd)
        print(kernel)     ---- this line prints information of all segments
    
    where:
        center - SPKID of central body (0 for SSB, 10 for Sun, etc.)
        target - SPKID of target body
        jd - time for computation (Julian date)
        position - a numpy array (x, y, z)
        velocity - a numpy array (xd, yd, zd)
        