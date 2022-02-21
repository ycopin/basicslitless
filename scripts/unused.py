def intensity2mag(I, refmag, refI):
    """Converts an intensity to a magnitude.

    Parameters
    ----------
    I : float
       intensity to convert.
    refmag : float
       reference magnitude.
    refI : float
       corresponding intensity

    Returns
    -------
    A magnitude corresponding to I
    """
    
    return refmag - 2.5*np.log10(I/refI)


def mag2intensity(mag, refI, refmag):
    """Converts a magnitude to an intensity.

    Parameters
    ----------
    mag : float
       magnitude to convert.
    refI : float
       reference intensity
    refmag : float
       corresponding magnitude.

    Returns
    -------
    An intensity corresponding to mag
    """
    
    return refI * 10**((refmag-mag)/2.5)


