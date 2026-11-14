#!/usr/bin/env python3

from decimal import Decimal, getcontext
from math import log2, log10, comb, e, factorial
from collections.abc import Hashable
import timeit


#### Generic caching layer to save time
class memoized(object):
  def __init__(self,func):
    self.func = func
    self.cache = {}
    self.__name__ = 'memoized:' + func.__name__
  def __call__(self,*args):
    if not isinstance(args, Hashable):
      return self.func(*args)
    if not args in self.cache:
      self.cache[args] = self.func(*args)
    return self.cache[args]


#### SPHINCS+ analysis

# Pr[FORS forgery given that exactly r sigs hit the leaf] = (1-(1-1/(2^a))^r)^k
@memoized
def forgeryprob(a, r, k):
  if k == 1: return 1 - (1 - 1 / Decimal(2**a))**r
  return forgeryprob(a, r, 1) * forgeryprob(a, r, k-1)

@memoized
def sigmalimit(tsec):
  return 1 / Decimal(2**tsec)

@memoized
def donelimit(tsec):
  return 1 - sigmalimit(tsec) / Decimal(2**20)

@memoized
def termlimit(tsec):
  return sigmalimit(tsec) / Decimal(2**20)

@memoized
def poisson_pmf(rate, n_events):
  return rate**n_events * Decimal(e)**(-rate) / factorial(n_events)

getcontext().prec = int(log10(2**768)) + 1  # ~232 digits: headroom for catastrophic cancellation

# Original SPHINCS+ paper's formula, improved by building the
# binomial pmf incrementally:  pmf_p(r) = pmf_p(r-1) * (sigs-r+1)/r * p/(1-p).
# This avoids recomputing comb(sigs, r) from scratch every step.
def compute_sigma_BHK19(leaves, a, k, sigs, tsec):
  sigma = Decimal(0)
  p = Decimal(1) / Decimal(leaves)
  q = Decimal(1) - p
  t = q ** sigs # pmf(0) = (1 - 1/leaves)**sigs
  done = t
  for r in range(1, sigs + 1):
    t = t * (Decimal(sigs) - r + 1) / Decimal(r) * (p / q)   # pmf(r)
    term = t * forgeryprob(a, r, k)
    sigma += term
    done += t
    if done >= donelimit(tsec) and term < termlimit(tsec):
      break
  return (sigma, r)

# Fluhrer & Dang's Poisson approximation.
def compute_sigma_FD24(leaves, a, k, sigs, tsec):
  sigma = Decimal(0)
  rate = Decimal(sigs) / Decimal(leaves)
  done = Decimal(poisson_pmf(rate, 0))
  for r in range(1, sigs + 1):
    t = poisson_pmf(rate, r)
    term = t * forgeryprob(a, r, k)
    sigma += term
    done += t
    if done >= donelimit(tsec) and term < termlimit(tsec):
      break
  return (sigma, r)

# Our contribution
def compute_sigma_ours(leaves, a, k, sigs, tsec):
  v = Decimal(1) - Decimal(1) / Decimal(2)**a
  t = lambda i: v**i
  sigma = Decimal(0)
  for i in range(k + 1):
    pgf = (Decimal(1) - (Decimal(1) - t(i)) / Decimal(leaves))**sigs
    sigma += comb(k, i) * ((-1)**i) * pgf
  return (sigma, k + 1)

# h, a, k, sigs, tsec
PARAMETER_SETS = {
  # FIPS-205
  'slh_dsa_128s':   [63, 12, 14, 2**64, 128],
  'slh_dsa_128f':   [66, 6,  33, 2**64, 128],
  'slh_dsa_192s':   [63, 14, 17, 2**64, 192],
  'slh_dsa_192f':   [66, 8,  33, 2**64, 192],
  'slh_dsa_256s':   [64, 14, 22, 2**64, 256],
  'slh_dsa_256f':   [68, 9,  35, 2**64, 256],

  # SP-800-230
  'slh_dsa_128_24': [22, 24, 6,  2**24, 128],
  'slh_dsa_192_24': [21, 25, 9,  2**24, 192],
  'slh_dsa_256_24': [21, 25, 12, 2**24, 256],

  # SHRINCS
  'shrincs': [45, 13, 10, 2**40, 128],

  # RLS
  'rls128gp1': [45, 23, 5, 2**30, 128],
  'rls128gp2': [45, 19, 6, 2**30, 128],
  'rls128gp3': [42, 24, 5, 2**30, 128],

  # FD24 divergence cases: small query count (low per-leaf hit prob) makes the
  # Poisson approximation diverge from the exact binomial.  Sorted by severity.
  'fd_diverge_1': [10,  9, 33, 4, 256],
  'fd_diverge_2': [10,  6, 33, 4, 256],
  'fd_diverge_3': [10,  9, 33, 8, 256],
  'fd_diverge_4': [10,  6, 33, 8, 256],
  'fd_diverge_5': [14,  6, 33, 4, 256],
  'fd_diverge_6': [10, 15, 20, 4, 256],

  # Performance comparison. BHK+19 and SD25 approximations converge very slowly on
  # this parameter set, requiring thousands of terms summed.
  'perf_stress': [52, 20, 16, 2**64, 128],
}

class _TimedCall:
  """Wraps a zero-arg callable so it can be timed by timeit while capturing
  its single return value (the wrapped function runs exactly once)."""
  def __init__(self, fn):
    self.fn = fn
    self.result = None
  def __call__(self):
    self.result = self.fn()


if __name__ == "__main__":
  for name, [h, a, k, sigs, tsec] in PARAMETER_SETS.items():
    # Time each computation with a single iteration (timeit), capturing its
    # result so the function is only run once per parameter set.
    our_call = _TimedCall(lambda: compute_sigma_ours(2**h, a, k, sigs, tsec))
    bhk_call = _TimedCall(lambda: compute_sigma_BHK19(2**h, a, k, sigs, tsec))
    fd_call = _TimedCall(lambda: compute_sigma_FD24(2**h, a, k, sigs, tsec))
    our_t = timeit.timeit(our_call, number=1)
    bhk_t = timeit.timeit(bhk_call, number=1)
    fd_t = timeit.timeit(fd_call, number=1)

    our_sigma, our_loops = our_call.result
    bhk_sigma, bhk_loops = bhk_call.result
    fd_sigma, fd_loops = fd_call.result

    bhk_err_pct = 100 * abs(our_sigma - bhk_sigma) / our_sigma
    fd_err_pct = 100 * abs(our_sigma - fd_sigma) / our_sigma

    print(name)
    print(f'  Ours (exact):        '
          f'sigma = 2**{log2(our_sigma):.19f} '
          f'({our_loops} loops, {our_t:.6f}s)')
    print(f'  BHK19 approximation: '
          f'sigma = 2**{log2(bhk_sigma):.19f} '
          f'({bhk_loops} loops, {bhk_err_pct:.8f}% error, {bhk_t:.6f}s)')
    print(f'  FD24 approximation:  '
          f'sigma = 2**{log2(fd_sigma):.19f} '
          f'({fd_loops} loops, {fd_err_pct:.8f}% error, {fd_t:.6f}s)')
