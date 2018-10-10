# -*- coding: utf-8 -*-
"""A supporting module for jplephem to handle data type 1 (Version 1.0.0)

This module computes position and velocity of a celestial small body, from a 
NASA SPICE SPK ephemeris kernel file of data type 1 (Modified Difference 
Arrays).
http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

You can get SPK files for many solar system small bodies from HORIZONS 
system of NASA/JPL.  See https://ssd.jpl.nasa.gov/?horizons

This module reads SPK files of data type 1, one of the types of binary SPK 
file.  

At the point of Oct. 2018, HORIZONS system provides files of type 21 for 
binary SPK files by default.  You can get type 1 binary SPK file for celestial 
small bodies through TELNET interface by answering back '1' for 
'SPK file format'.

Module required:
    jplephem (version 2.6 or later)
    numpy

Usage:
    from spktype01 import SPKType01
    kernel = SPKType01.open('path')
    position, velocity = kernel.compute_type01(center, target, jd)
    
    where:
        center - SPKID of central body (0 for SSB, 10 for Sun, etc.)
        target - SPKID of target body
        jd - time for computation (Julian date)

Author: Shushi Uetsuki (whiskie14142)
This module has been developed based on jplephem.spk and FORTRAN source 
of the SPICE Toolkit of NASA/JPL/NAIF.
jplephem : https://pypi.python.org/pypi/jplephem
SPICE Toolkit : http://naif.jpl.nasa.gov/naif/toolkit.html
"""

from numpy import array, zeros, reshape
from jplephem.daf import DAF
from jplephem.names import target_names

T0 = 2451545.0
S_PER_DAY = 86400.0

def jd(seconds):
    """Convert a number of seconds since J2000 to a Julian Date.
    """
    return T0 + seconds / S_PER_DAY

class SPKType01(object):
    """Class for SPK kernel to handle data type 1 (Modified Difference Arrays)
    """
    def __init__(self, daf):
        self.daf = daf
        self.segments = [Segment(self.daf, *t) for t in self.daf.summaries()]
        ssec = lambda s : s.start_second
        self.segments.sort(key=ssec)
        
        # initialize arrays for spke01
        self.G = zeros(15)
        self.REFPOS = zeros(3)
        self.REFVEL = zeros(3)
        self.KQ = array([0, 0, 0])
        self.FC = zeros(15)
        self.FC[0] = 1.0
        self.WC = zeros(13)
        self.W = zeros(17)
        
        # initialize for compute_type01
        self.mda_record_exist = False
        self.current_segment_exist = False
        

    @classmethod
    def open(cls, path):
        """Open the file at `path` and return an SPK instance.
        """
        return cls(DAF(open(path, 'rb')))

    def close(self):
        """Close this SPK file."""
        self.daf.file.close()

    def __str__(self):
        daf = self.daf
        d = lambda b: b.decode('latin-1')
        lines = (str(segment) for segment in self.segments)
        return 'File type {0} and format {1} with {2} segments:\n{3}'.format(
            d(daf.locidw), d(daf.locfmt), len(self.segments), '\n'.join(lines))
    
    def comments(self):
        return self.daf.comments()

    def compute_type01(self, center, target, jd1, jd2=0.0):
        """Compute position and velocity of target from SPK data (data type 1).
        Inputs:
            center - SPKID of the central body (0 for SSB, 10 for Sun, etc)
            target - SPKID of the target
            jd1, jd2 - Julian date of epoch for computation.  (jd1 + jd2) will 
                be used for computation.  If you want precise definition of 
                epoch, jd1 should be an integer or a half integer, and jd2 should
                be a relatively small floating point number.
        Returns:
            Position (X, Y, Z) and velocity (XD, YD, ZD) of the target at 
            epoch.  Position and velocity are provided as Numpy arrays 
            respectively.
        """
        eval_sec = (jd1 - T0)
        eval_sec = (eval_sec + jd2) * S_PER_DAY
        
        if self.mda_record_exist:
            if eval_sec >= self.mda_lb and eval_sec < self.mda_ub:
                result = self.spke01(eval_sec, self.mda_record)
                return result[0:3], result[3:]
        
        self.mda_record, self.mda_lb, self.mda_ub = self.get_MDA_record(eval_sec, target, center)
        self.mda_record_exists = True
        
        result = self.spke01(eval_sec, self.mda_record)
        return result[0:3], result[3:]
                
    def get_MDA_record(self, eval_sec, target, center):
        """Return a MDA record for defined epoch.
        Inputs:
            eval_sec - epoch for computation, seconds from J2000
            target - body ID of the target
            center - body ID of coordinate center
        Returns:
            MDA record - a Numpy array of 71 floating point numbers
        Exception:
            ValueError will be raised when:
                eval_sed is outside of SPK data
                target and center are not in SPK data
                invalid data type of SPK data
        """
        
        # chech last segment can be used
        if self.current_segment_exist:
            if eval_sec >= self.current_segment.start_second    \
                and eval_sec < self.current_segment.end_second  \
                and target == self.current_segment.target            \
                and center == self.current_segment.center:
                
                return self.current_segment.get_MDA_record(eval_sec)

        # select segments with matched 'target' and 'center'
        matched = []
        for segment in self.segments:
            if segment.target == target and segment.center == center:
                matched.append(segment)
        if len(matched) == 0:
            raise ValueError('Invalid Target and/or Center')
        if eval_sec < matched[0].start_second or eval_sec >= matched[-1].end_second:
            raise ValueError('Invalid Time to evaluate')
        
        # selet a segment based on eval_sec
        found = False
        for segment in matched:
            if eval_sec < segment.end_second:
                found = True
                self.current_segment = segment
                break
        if not found:
            self.current_segment = matched[-1]
        self.current_segment_exist = True
        
        # get the MDA record from selected segment
        if self.current_segment.data_type != 1:
            raise ValueError('Invalid data. Data Type must be 1.')
        
        return self.current_segment.get_MDA_record(eval_sec)
        
    def spke01(self, ET, RECORD):
        """Compute position and velocity from a Modified Difference Array record
        
        Inputs:
            ET: Epoch time to evaluate position and velocity (seconds since J2000)
            RECORD: A record of Modified Difference Array
        Returns: STATE
            STATE: A numpy array which contains position and velocity
        """
        
        # This method has been translated from SPKE01 of SPICE Toolkit and
        # modified by Shushi Uetsuki.
        #
        # SPICE Toolkit for FORTRAN : http://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html
        # SPK Required Reading : http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
        #
        # Original FORTRAN code uses 'SAVE' directive, and it means all variable
        # should be saved for next call.  So i decided to make almost all 
        # variables to be instance variable.  Some of them are initialized in 
        # __init__ method.
        
        STATE = zeros(6)

#     Variable  I/O  Description
#     --------  ---  --------------------------------------------------
#     ET         I   Target epoch.
#     RECORD     I   Data record.
#     STATE      O   State (position and velocity).
#
#$ Detailed_Input
#
#     ET          is a target epoch, at which a state vector is to
#                 be computed.
#
#     RECORD      is a data record which, when evaluated at epoch ET,
#                 will give the state (position and velocity) of some
#                 body, relative to some center, in some inertial
#                 reference frame.
#
#$ Detailed_Output
#
#     STATE       is the state. Units are km and km/sec.
#
#$ Parameters
#
#     None.
#
#$ Exceptions
#
#     None.
#
#$ Files
#
#     None.
#
#$ Particulars
#
#     The exact format and structure of type 1 (difference lines)
#     segments are described in the SPK Required Reading file.
#
#     Difference lines (DL's) are generated by JPL navigation
#     system programs P and PV. Each data record is equivalent
#     to the (slightly rearranged) 'P' portion of a NAVIO PV file
#     data record.
#
#     SPKE01 is a specialized version of Fred Krogh's subroutine DAINT.
#     Only the calling sequence has been changed.
#
#     Because the original version was undocumented, only Fred 
#     knows how this really works.
#
#$ Examples
#
#     None.
#
#$ Restrictions
#
#     Unknown.
#
#$ Literature_References
#
#     NAIF Document 168.0, "S- and P- Kernel (SPK) Specification and
#     User's Guide"
#
#$ Author_and_Institution
#
#     F.T. Krogh      (JPL)
#     I.M. Underwood  (JPL)
#
#$ Version
#
#-    SPICELIB Version 1.1.0, 14-FEB-1997 (WLT) 
#     
#        The goto's were removed and loop and if structures 
#        revealed.  We still don't know exactly what's going 
#        on, but at least the bones of this routine have been 
#        cleaned off and are ready for assembly. (WLT)
#
#-    SPICELIB Version 1.0.4, 30-OCT-1996 (WLT)
#
#        Removed redundant SAVE statements from the declaration
#        section.  Thanks to Steve Schlaifer for finding this
#        error.
#
#-    SPICELIB Version 1.0.3, 10-MAR-1992 (WLT)
#
#        Comment section for permuted index source lines was added
#        following the header.
#
#-    SPICELIB Version 1.0.2, 23-AUG-1991 (HAN)
#
#        SPK01 was removed from the Required_Reading section of the
#        header. The information in the SPK01 Required Reading file
#        is now part of the SPK Required Reading file.
#
#-    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN)
#
#        Literature references added to the header.
#
#-    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) (FTK)
#
#-&
# 
#$ Index_Entries
#
#     evaluate type_1 spk segment
#
#-&


        
#
#     Unpack the contents of the MDA array.
#
#        Name    Dimension  Description
#        ------  ---------  -------------------------------
#        TL              1  Final epoch of record
#        G              15  Stepsize function vector
#        REFPOS          3  Reference position vector
#        REFVEL          3  Reference velocity vector
#        DT         15,NTE  Modified divided difference arrays
#        KQMAX1          1  Maximum integration order plus 1
#        KQ            NTE  Integration order array
#
#     For our purposes, NTE is always 3.
#
        self.TL = RECORD[0]
        self.G = RECORD[1:1+15]
#     
#     Collect the reference position and velocity.
#     
        self.REFPOS[0] = RECORD[16]
        self.REFVEL[0] = RECORD[17]
        
        self.REFPOS[1] = RECORD[18]
        self.REFVEL[1] = RECORD[19]
        
        self.REFPOS[2] = RECORD[20]
        self.REFVEL[2] = RECORD[21]
        
        self.DT = reshape(RECORD[22:22+45], (15, 3), order='F')
        
        self.KQMAX1 = int(RECORD[67])
        self.KQ[0] = int(RECORD[68])
        self.KQ[1] = int(RECORD[69])
        self.KQ[2] = int(RECORD[70])
#     
#     Next we set up for the computation of the various differences
#     
        self.DELTA = ET - self.TL
        self.TP = self.DELTA
        self.MQ2 = self.KQMAX1 - 2
        self.KS = self.KQMAX1 - 1
#
#     This is clearly collecting some kind of coefficients.  
#     The problem is that we have no idea what they are...
#     
#     The G coefficients are supposed to be some kind of step size 
#     vector. 
#     
#     TP starts out as the delta t between the request time 
#     and the time for which we last had a state in the MDL file. 
#     We then change it from DELTA  by the components of the stepsize 
#     vector G.  
#
        for J in range(1, self.MQ2 + 1):
            self.FC[J] = self.TP / self.G[J-1]
            self.WC[J-1] = self.DELTA / self.G[J-1]
            self.TP = self.DELTA + self.G[J-1]
#
#     Collect KQMAX1 reciprocals. 
#   
        for J in range(1, self.KQMAX1 + 1):
            self.W[J-1] = 1.0 / float(J)
#
#     Compute the W(K) terms needed for the position interpolation
#     (Note,  it is assumed throughout this routine that KS, which 
#     starts out as KQMAX1-1 (the ``maximum integration'') 
#     is at least 2.
#
        self.JX = 0
        self.KS1 = self.KS - 1
        
        while self.KS >= 2:
            
            self.JX = self.JX + 1
            
            for J in range(1, self.JX + 1):
                self.W[J+self.KS-1] = self.FC[J] * self.W[J+self.KS1-1] - self.WC[J-1] * self.W[J+self.KS-1]
            
            self.KS = self.KS1
            self.KS1 = self.KS1 - 1
#
#     Perform position interpolation: (Note that KS = 1 right now.
#     We don't know much more than that.)
#
        for I in range(1, 3 + 1):
            
            self.KQQ = self.KQ[I-1]
            self.SUM = 0.0
            
            for J in range(self.KQQ, 0, -1):
                self.SUM = self.SUM + self.DT[J-1, I-1] * self.W[J+self.KS-1]
            
            STATE[I-1] = self.REFPOS[I-1] + self.DELTA * (self.REFVEL[I-1] + self.DELTA * self.SUM)
#
#     Again we need to compute the W(K) coefficients that are 
#     going to be used in the velocity interpolation. 
#     (Note, at this point, KS = 1, KS1 = 0.)
#      
        for J in range(1, self.JX + 1):
            self.W[J+self.KS-1] = self.FC[J] * self.W[J+self.KS1-1] - self.WC[J-1] * self.W[J+self.KS-1]
        
        self.KS = self.KS - 1
        
#
#     Perform velocity interpolation:
#
        for I in range(1, 3 + 1):
            self.KQQ = self.KQ[I-1]
            self.SUM = 0.0
            
            for J in range(self.KQQ, 0, -1):
                self.SUM = self.SUM + self.DT[J-1, I-1] * self.W[J+self.KS-1]
            
            STATE[I+3-1] = self.REFVEL[I-1] + self.DELTA * self.SUM
        
#
#     That's all folks.  We don't know why we did anything, but 
#     at least we can tell structurally what we did.
#      
        
        return STATE
        
        

class Segment(object):
    """A single segment of a SPK file.

    There are several items of information about each segment that are
    loaded from the underlying SPK file, and made available as object
    attributes:

    segment.source - official ephemeris name, like 'DE-0430LE-0430'
    segment.start_second - initial epoch, as seconds from J2000
    segment.end_second - final epoch, as seconds from J2000
    segment.start_jd - start_second, converted to a Julian Date
    segment.end_jd - end_second, converted to a Julian Date
    segment.center - integer center identifier
    segment.target - integer target identifier
    segment.frame - integer frame identifier
    segment.data_type - integer data type identifier
    segment.start_i - index where segment starts
    segment.end_i - index where segment ends
    """
    def __init__(self, daf, source, descriptor):
        self.daf = daf
        self.source = source
        (self.start_second, self.end_second, self.target, self.center,
         self.frame, self.data_type, self.start_i, self.end_i) = descriptor
        self.start_jd = jd(self.start_second)
        self.end_jd = jd(self.end_second)

    def __str__(self):
        return self.describe(verbose=False)

    def describe(self, verbose=True):
        """Return a textual description of the segment.
        """
        center = titlecase(target_names.get(self.center, 'Unknown center'))
        target = titlecase(target_names.get(self.target, 'Unknown target'))
        text = ('{0.start_jd:.2f}..{0.end_jd:.2f}  {1} ({0.center})'
                ' -> {2} ({0.target})'
                ' data_type={0.data_type}'.format(self, center, target))
        if verbose:
            text += ('\n  frame={0.frame} data_type={0.data_type} source={1}'
                     .format(self, self.source.decode('ascii')))
        return text
        
    def get_MDA_record(self, time_sec):
        """Return a Modified Difference Array(MDA) record for the time to 
        evaluate with its effective time boundaries (lower and upper).
        Inputs:
            time_sec - epoch for computation, seconds from J2000
        Returns: mda_record, lower_boundary, upper_boundary
            mda_record: A Modified Difference Array record
            lower_boundary: lower boundary of the record, seconds since J2000
            upper_boundary: upper boundary of the record, seconds since J2000
        """

        # Number of records in this segment
        entry_count = int(self.daf.map_array(self.end_i, self.end_i))
        
        # Number of entries in epoch directory 
        epoch_dir_count = entry_count // 100
        
        # serch target epoch in epoch directory to narrow serching aria
        if epoch_dir_count >= 1:
            epoch_dir = self.daf.map_array(self.end_i - epoch_dir_count,
                                            self.end_i - 1)
            found = False
            for i in range(1, epoch_dir_count + 1):
                if epoch_dir[i-1] > time_sec:
                    found = True
                    break
            if found:
                serch_last_index = i * 100
                serch_start_index = (i - 1) * 100 + 1
            else:
                serch_last_index = entry_count
                serch_start_index = epoch_dir_count * 100 + 1
        else:
            serch_last_index = entry_count
            serch_start_index = 1

        # epoch_table contains epochs for all records in this segment        
        epoch_table = self.daf.map_array(self.start_i + (entry_count * 71),
                                       self.start_i + (entry_count * 71) + entry_count - 1)

        # serch target epoch in epoch_table
        found = False
        for i in range(serch_start_index, serch_last_index + 1):
            if epoch_table[i-1] > time_sec:
                found = True
                break
        if not found:
            i = serch_last_index
        record_index = i
        upper_boundary = epoch_table[i-1]
        if i != 1:
            lower_boundary = epoch_table[i-2]
        else:
            lower_boundary = self.start_second
        
        mda_record = self.daf.map_array(self.start_i + ((record_index - 1) * 71),
                                        self.start_i + (record_index * 71) - 1)

        # mda_record : one record of MDA
        # lower_boundary : lower boundary of epoch in this MDA record
        # upper_boundary : upper boundary of epoch in this MDA record
        return mda_record, lower_boundary, upper_boundary

def titlecase(name):
    """Title-case target `name` if it looks safe to do so.
    """
    return name if name.startswith(('1', 'C/', 'DSS-')) else name.title()







