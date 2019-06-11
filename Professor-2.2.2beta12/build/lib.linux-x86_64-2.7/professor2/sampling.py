# -*- python -*-



# http://stackoverflow.com/questions/3025162/statistics-combinations-in-python
def numCombs(n, k):
    """
    n choose k algorithm
    """
    from operator import mul
    from fractions import Fraction
    return int( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )

def xrandomUniqueCombinations(items, nchoose, howmany=None):
    """ Generator-like function for n choose k items """
    seencombs = []
    # Max number safeguard against infinite loops
    maxnum = numCombs(len(items), nchoose)
    import random
    if howmany is None or howmany > maxnum:
        print "Only %i possible combinations"%maxnum
        howmany = maxnum
    while len(seencombs) < howmany:
        temp = random.sample(items, nchoose)
        temp.sort()
        if not sorted(temp) in seencombs:
            seencombs.append(temp)
            yield temp


## Define a sampler type
# This is obsolete now
class Sampler(object):
    # @deprecated
    def __init__(self, low, high, bias=None):
        self.low = float(low)
        self.high = float(high)
        self.f, self.invf = None, None
        if bias:
            ## Import clever machinery
            try:
                import sympy as sp
                from sympy.abc import x, y
                import numpy as np
            except ImportError:
                print "Bias functions require SymPy and NumPy to be installed... exiting"
                exit(1) #< TODO: don't exit from inside a lib function...
            ## Make transformation and its inverse
            try:
                #print bias
                f_expr = sp.sympify(bias)
            except sp.SympifyError, e:
                print "Bias function could not be parsed by SymPy:"
                print e
                exit(1) #< TODO: don't exit from inside a lib function...
            try:
                finv_exprs = sp.solve(sp.Eq(y, f_expr), x)
                finv_expr = finv_exprs[0]
                #print f_expr, finv_exprs
                self.f = sp.lambdify(x, f_expr, "numpy")
                self.finv = sp.lambdify(y, finv_expr, "numpy")
                self.lowf, self.highf = self.f(self.low), self.f(self.high)
            except Exception, e:
                print "Bias function could not be used/inverted by SymPy:"
                print e
                exit(1) #< TODO: don't exit from inside a lib function...

    def shoot(self):
        import random
        if not self.f:
            ## Just uniform sampling between low..high
            val = random.uniform(self.low, self.high)
        else:
            ## Uniform sample in transformed space, and transform the result back
            valf = random.uniform(self.lowf, self.highf)
            val = self.finv(valf)
        return val

    def __call__(self):
        return self.shoot()

    def __repr__(self):
        return "<%s with x in %f ... %f>" % (self.__class__.__name__, self.low, self.high)



class Sobol(object):

    def __init__(self, dim, initialSeed=0):
        self._seed=initialSeed
        self._dim=dim

    def shoot(self):
        try:
            import sobol
        except ImportError:
            print "sobol not available, try pip install sobol"
            exit(1)
        p, newseed = sobol.sobol_seq.i4_sobol(self._dim, self._seed)
        self._seed=newseed
        return p

    def __call__(self):
        return self.shoot()

class RandomU(object):

    def __init__(self, dim):
        self._dim=dim

    def shoot(self):
        from numpy import random
        return random.uniform(0,1,(self._dim,1)).flatten()

    def __call__(self):
        return self.shoot()

class NDSampler(object):

    def __init__(self, ranges, biases=None, sobol=False, seed=None):
        try:
            from collections import OrderedDict
        except:
            from ordereddict import OrderedDict
        self._ranges = ranges
        self._dim = len(self._ranges)
        self._sobol=sobol
        self._biases=biases
        self.f, self.invf = [None for i in xrange(self._dim)], [None for i in xrange(self._dim)]

        if seed is not None:
            print "Warning: Seed setting currently not implemented"

        if biases is not None:
            self.setBiases()

        if sobol:
            self._generator=Sobol(self._dim)
        else:
            self._generator=RandomU(self._dim)

    def setBiases(self):
        try:
            import sympy as sp
            from sympy.abc import x, y
            import numpy as np
        except ImportError:
            print "Bias functions require SymPy and NumPy to be installed"
            print "Try pip install sympy numpy"
            exit(1) #< TODO: don't exit from inside a lib function...
            # self.lowf, self.highf = self.f(self.low), self.f(self.high)
        for num, b in enumerate(self._biases):
            if b is not None:
                ## Make transformation and its inverse
                try:
                    f_expr = sp.sympify(b)
                except sp.SympifyError, e:
                    print "Bias function could not be parsed by SymPy:"
                    print e
                    exit(1) #< TODO: don't exit from inside a lib function...
                try:
                    finv_exprs = sp.solve(sp.Eq(y, f_expr), x)
                    finv_expr = finv_exprs[0]
                    #print f_expr, finv_exprs
                    self.f[num] = sp.lambdify(x, f_expr, "numpy")
                    self.invf[num] = sp.lambdify(y, finv_expr, "numpy")
                except Exception, e:
                    print "Bias function could not be used/inverted by SymPy:"
                    print e
                    exit(1) #< TODO: don't exit from inside a lib function...

    def scale(self, Praw):
        P=[]
        for num, p in enumerate(Praw):
            a, b = self._ranges[num]
            if self.f[num] is not None:
                a=self.f[num](a)
                b=self.f[num](b)
            pscaled = a + p*(b-a)
            if self.invf[num] is not None:
                pscaled = self.invf[num](pscaled)
            P.append(pscaled)
        return P


    def __call__(self):
        P_raw = self._generator.shoot()
        return self.scale(P_raw)

    def __repr__(self):
        s="<%s,  %i-D %s>"%(self.__class__.__name__, self._dim, "Sobol" if self._sobol else "Uniform")
        return s


## Test biased sampler machinery if run as main
if __name__ == "__main__":
    s = Sampler(1, 10, "exp(x)")
    import yoda
    h = yoda.Histo1D(20, 0, 10)
    for _ in xrange(10000):
        h.fill( s() )
    yoda.plot(h, "foo.pdf")

    import pylab
    s=Sobol(2)
    PS = [s() for _ in xrange(50)]
    X=[p[0] for p in PS]
    Y=[p[1] for p in PS]
    pylab.clf()
    pylab.plot(X,Y, "bo", label="Sobol")

    s=RandomU(2)
    PU = [s() for _ in xrange(50)]
    X=[p[0] for p in PU]
    Y=[p[1] for p in PU]
    pylab.plot(X, Y, "rx", label="Random uniform")
    pylab.legend()
    pylab.savefig("foobol.pdf")

    N = NDSampler([[1,10]], ["exp(x)"], True)
    h = yoda.Histo1D(20, 0, 10)
    for _ in xrange(10000):
        h.fill( N()[0] )

    pylab.clf()
    yoda.plot(h, "foobolN.pdf")



    pylab.clf()
    NS = NDSampler([[1,10], [2,3]], ["exp(x)", "10**x"], True)
    print NS
    P2S = [NS() for _ in xrange(50)]
    X=[p[0] for p in P2S]
    Y=[p[1] for p in P2S]
    pylab.plot(X, Y, "bo", label="Sobol")

    NU = NDSampler([[1,10], [2,3]], ["exp(x)", "10**x"], False)
    print NU
    P2U = [NU() for _ in xrange(50)]
    X=[p[0] for p in P2U]
    Y=[p[1] for p in P2U]
    pylab.plot(X, Y, "rx", label="Random uniform")

    pylab.legend()
    pylab.savefig("foosobol2_wbias.pdf")

