# -*- python -*-

def mk_fitfunc(fname, pnames, globname, extraargs=[]):
    """
    Dynamically make a fit function for the given param names, to be passed to Minuit.

    Return a string definition of the function, to be exec'd, and the list of
    generated internal arg names corresponding to pnames.
    """
    fargs = ["A%03i" % i for i in xrange(len(pnames))]
    funcdef = "def {gname}({fargs}):\n    x= lambda {fargs}: {fname}([{fargs}], {extrargs})".format(gname=globname, fargs=", ".join(fargs), fname=fname, extrargs=", ".join(extraargs))
    funcdef+="\n    return x({fargs})".format(fargs=", ".join(fargs))
    return funcdef


def prepareBins(dataHistos, ipolHistos, maxErrDict, matchers, doFilter=False):
    ## Try to read run histos and extract maximum errors

    ## Weight file parsing
    # ?

    ## Find things available in both ipol and ref data, and in the weight file if there is one
    available = []
    for ihn in sorted(ipolHistos.keys()):
        ## Set default bin weights
        for ib in ipolHistos[ihn].bins:
            ib.w = 1.0
        ## Find user-specified bin weights if there was a weight file
        if matchers is not None:
            ## Find matches
            pathmatch_matchers = [(m,wstr) for m,wstr in matchers.iteritems() if m.match_path(ihn)]
            ## Ditch histos not listed in the weight file
            if not pathmatch_matchers:
                del ipolHistos[ihn]
                continue
            ## Attach fit weights to the ibins, setting to zero if there's no position match
            for ib in ipolHistos[ihn].bins:
                posmatch_matchers = [(m,wstr) for (m,wstr) in pathmatch_matchers if m.match_pos(ib)]
                ib.w = float(posmatch_matchers[-1][1]) if posmatch_matchers else 0 #< NB. using last match
        for rhn in dataHistos.keys():
            if ihn == rhn or rhn == "/REF/"+ihn: #< TODO: short for rhn = "/REF/"+ihn ?
                # TODO: we should eliminate this potential mismatch of ref and MC hnames
                available.append([ihn,rhn])
                break #< TODO: ok?
            # else:
                # print "Could not find %s"%ihn


    ## Prepare lists of ibins and dbins
    IBINS, DBINS, MAXERRS, FILTERED, INCONSISTENT = [], [], [], [], []
    BINDICES = {} # Allows for more helpful error messages in case of prof.StatError
    for a in available:
        # TODO: print out the available observables
        if len(ipolHistos[a[0]].bins) != len(dataHistos[a[1]].bins):
            INCONSISTENT.append(a[0])
            del ipolHistos[a[0]]
            del dataHistos[a[1]]
        else:
            BINDICES[a[0]] = []#range(len(IBINS),  len(IBINS) +     len(ipolHistos[a[0]])) # This is for debugging
            for nb in xrange(len(ipolHistos[a[0]].bins)):
                if doFilter and dataHistos[a[1]].bins[nb].err ==0:
                    FILTERED.append(1)
                    continue
                if ipolHistos[a[0]].bins[nb].w >0:
                    IBINS.append(ipolHistos[a[0]].bins[nb])
                    DBINS.append(dataHistos[a[1]].bins[nb])
                    BINDICES[a[0]].append(len(IBINS))
            if maxErrDict:
                MAXERRS.extend(MAXERRS[a[0]])

    if not MAXERRS:
        MAXERRS = None
    if INCONSISTENT:
        print "Inconsistencies discovered between data bins and parametrised bins:"
        for h in INCONSISTENT:
            print "  Removing histogram", h
    if FILTERED:
        print "Filtered %i bins due to zero data error" % len(FILTERED)



    ## Sanity checks
    assert len(IBINS) == len(DBINS)
    if not IBINS:
        print "No bins... exiting"
        import sys
        sys.exit(1)
    assert MAXERRS is None or len(IBINS) == len(MAXERRS)
    return DBINS, IBINS, MAXERRS


def simpleGoF(params, dbins, ibins, maxerrs, extraErr=1.):
    """
    Very straightforward goodness-of-fit measure
    """
    chi2 = 0.0
    for num, ibin in enumerate(ibins):
        ## Weight is attached to the ipol bin (default set to 1.0 above)
        w = ibin.w
        if w == 0:
            continue
        ## Get ipol & ref bin values and compute their difference
        ival = ibin.val(params)
        dval = dbins[num].val
        diff = dval - ival
        ## Data error
        err2 = dbins[num].err**2
        ## Plus interpolation error added in quadrature
        maxierr = maxerrs[ibin] if maxerrs else None
        err2 += ibin.err(params, emax=maxierr)**2
        # TODO: compute asymm error for appropriate deviation direction cf. sum([e**2 for e in ibin.ierrs])
        if not err2:
            continue
        # TODO: should we square w too, so it penalised deviations _linearly_?
        chi2 += w * diff**2 / err2
    return chi2#/extraErr**2


# TODO make this work again
def simpleGoFGradient(params):
    """
    Very straightforward goodness-of-fit measure
    """
    dchi2 = [0 for x in params]
    for num, ibin in enumerate(IBINS):
        ## Weight is attached to the ipol bin (default set to 1.0 above)
        w = ibin.w
        if w == 0:
            continue
        ## Get ipol & ref bin values and compute their difference
        ival = ibin.val(params)
        dval = DBINS[num].val
        diff = dval - ival
        ## Data error
        err2 = DBINS[num].err**2
        ## Plus interpolation error added in quadrature
        maxierr = MAXERRS[ibin] if MAXERRS else None
        err2 += ibin.err(params, emax=maxierr)**2
        # TODO: compute asymm error for appropriate deviation direction cf. sum([e**2 for e in ibin.ierrs])
        if not err2:
            raise prof.StatError("Zero uncertainty on a bin being used in the fit -- cannot compute a reasonable GoF")
        # TODO: should we square w too, so it penalised deviations _linearly_?
        igrad = ibin.grad(params)
        for p in xrange(len(params)):
            dchi2[p] += 2 * w * diff * igrad[p] / err2
    N=sum(dchi2)
    return [x/N for x in dchi2]

def setupMinuitFitarg(pnames, pmins, pmaxs, LIMITFILE):
    ## Dictionary fitarg for iminuit
    farg=dict()

    ## Initial conditions --- use pos = center of hypercube, and step = range/10
    # TODO: Optionally make an initial brute force scan to choose the Minuit starting point, using prof.scangrid
    assert len(pmins) == len(pmaxs)

    pmids = [(pmins[i] + pmaxs[i])/2. for i in xrange(len(pmins))]
    pranges = [(pmaxs[i] - pmins[i]) for i in xrange(len(pmins))]

    # This sets the start point
    for i, aname in enumerate(pnames):
        farg[aname] = pmids[i]
        farg['error_%s'%aname] = pranges[i] / 10.

    ## Fix parameters, set limits (with pname translation)
    import professor2 as prof
    limits, fixed = prof.read_limitsandfixed(LIMITFILE)

    for i, pname in enumerate(pnames):
        if pname in limits.keys():
            farg['limit_%s'%pname] = limits[pname]
        if pname in fixed.keys():
            # if not opts.QUIET:
                # print "Fixing", pname, "= %f"%fixed[pnames[i]]
            farg[pname] = fixed[pnames[i]]
            farg['fix_%s'%pname] = True
    return farg


def readResult(fname):
    """
       Open results file, extract and return minimum point as OrderedDict
       and return raw list of all other lines for further processing.
    """
    RES=[]
    OTH=[]
    with open(fname) as f:
        for line in f:
            l=line.strip()
            if l.startswith("#"):
                OTH.append(l)
            else:
                temp=l.split()
                RES.append([temp[0], float(temp[1])])

    from collections import OrderedDict
    return OrderedDict(RES), OTH

def getParamCov(TXT):
    """
       Read the covariance matrix from the lines, return as numpy array
    """
    START = TXT.index("# Covariance matrix:") + 2
    dim = len(TXT[START].strip().split()) - 2
    END = START+dim
    COV_raw = TXT[START:END]
    from numpy import zeros
    COV_p = zeros((dim, dim))
    for i in xrange(dim):
        temp = map(float, COV_raw[i].split()[2:2+dim])
        for j in xrange(dim):
            COV_p[i][j] = temp[j]


    return COV_p

def getParamCorr(TXT):
    """
       Read the corrlation matrix from the lines, return as numpy array
    """
    START = TXT.index("# Correlation matrix:") + 2
    dim = len(TXT[START].strip().split()) - 2
    END = START+dim
    COV_raw = TXT[START:END]
    from numpy import zeros
    COV_p = zeros((dim, dim))
    for i in xrange(dim):
        temp = map(float, COV_raw[i].split()[2:2+dim])
        for j in xrange(dim):
            COV_p[i][j] = temp[j]


    return COV_p

def eigenDecomposition(mat):
    """
    Given a symmetric, real NxN matrix, M, an eigen decomposition is always
    possible, such that M can be written as M = T_transp * S * T (orthogonal
    transformation) with T_transp * T = 1_N and S being a diagonal matrix with
    the eigenvalues of M on the diagonal.

    Returns
    -------
    T_trans : numpy.matrix
    S : numpy.ndarray
        The real eigen values.
    T : numpy.matrix
    """
    from scipy import linalg
    from numpy import matrix
    import numpy

    A = matrix(mat)
    S, T_trans = linalg.eig(A)
    if numpy.iscomplex(S).any():
        raise ValueError("Given matrix `mat` has complex eigenvalues!")

    return matrix(T_trans), S.real, matrix(T_trans).transpose()

def mkEigentunes(COV, point, plus=True):
    T_trans, S, T = eigenDecomposition(COV)
    from numpy import sqrt, zeros, matrix
    rv = matrix(point.values())
    rv_trans = (T_trans * rv.transpose()).transpose()

    ret = []

    for num, c in enumerate(S):
        ev = zeros(len(S))
        sigma=sqrt(c)
        if plus:
            ev[num] =sigma
        else:
            ev[num] = -1* sigma
        ev_trans = rv_trans + ev
        etune_params_t = T * ev_trans.transpose()
        etune_params = etune_params_t.transpose().tolist()[0]
        ret.append([sigma, etune_params])

    return ret
