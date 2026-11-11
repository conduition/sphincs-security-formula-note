#!/usr/bin/env python3

from decimal import Decimal, getcontext
from math import log2, log10, comb, e, factorial


#### SPHINCS+ analysis

# Pr[exactly r sigs hit the leaf targeted by this forgery attempt]
def qhitprob(leaves, qs, r):
  p = 1 / Decimal(leaves)
  return comb(qs, r) * p**r * (1 - p)**(qs - r)

# Pr[FORS forgery given that exactly r sigs hit the leaf] = (1-(1-1/(2^a))^r)^k
def forgeryprob(a, r, k):
  if k == 1: return 1 - (1 - 1 / Decimal(2**a))**r
  return forgeryprob(a, r, 1) * forgeryprob(a, r, k-1)

def sigmalimit(tsec):
  return 1 / Decimal(2**tsec)

def donelimit(tsec):
  return 1 - sigmalimit(tsec) / Decimal(2**20)

def poisson_pmf(rate, n_events):
  return rate**n_events * Decimal(e)**(-rate) / factorial(n_events)

getcontext().prec = int(log10(2**768)) + 1  # ~232 digits: headroom for catastrophic cancellation


# Original SPHINCS+ paper's formula
def compute_sigma_BHK19(leaves, a, k, sigs, tsec):
  sigma = Decimal(0)
  done = qhitprob(leaves, sigs, 0)
  for r in range(1, sigs + 1):
    t = qhitprob(leaves, sigs, r)
    term = t * forgeryprob(a, r, k)
    sigma += term
    done += t
    if done >= donelimit(tsec) and term < sigmalimit(tsec) / Decimal(2**20):
      break
  return (sigma, r)

# Fluhrer & Dang's Poisson approximation.
def compute_sigma_FD25(leaves, a, k, sigs, tsec):
  sigma = Decimal(0)
  rate = Decimal(sigs) / Decimal(leaves)
  done = Decimal(poisson_pmf(rate, 0))
  for r in range(1, sigs + 1):
    t = poisson_pmf(rate, r)
    term = t * forgeryprob(a, r, k)
    sigma += term
    done += t
    if done >= donelimit(tsec) and term < sigmalimit(tsec) / Decimal(2**20):
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

  # FD25 divergence cases: small hypertree (large per-leaf hit prob) makes the
  # Poisson approximation diverge from the exact binomial.  Sorted by severity.
  'fd_diverge_1': [10,  9, 33, 4, 256],
  'fd_diverge_2': [10,  6, 33, 4, 256],
  'fd_diverge_3': [10,  9, 33, 8, 256],
  'fd_diverge_4': [10,  6, 33, 8, 256],
  'fd_diverge_5': [14,  6, 33, 4, 256],
  'fd_diverge_6': [10, 15, 20, 4, 256],
}

if __name__ == "__main__":
  for name, [h, a, k, sigs, tsec] in PARAMETER_SETS.items():
    our_sigma, our_loops = compute_sigma_ours(2**h, a, k, sigs, tsec)
    bhk_sigma, bhk_loops = compute_sigma_BHK19(2**h, a, k, sigs, tsec)
    fd_sigma, fd_loops = compute_sigma_FD25(2**h, a, k, sigs, tsec)

    bhk_err_pct = 100 * abs(our_sigma - bhk_sigma) / our_sigma
    fd_err_pct = 100 * abs(our_sigma - fd_sigma) / our_sigma

    print(name)
    print(f'  Ours (exact):        sigma = 2**{log2(our_sigma):.19f} ({our_loops} loops)')
    print(f'  BHK19 approximation: sigma = 2**{log2(bhk_sigma):.19f} ({bhk_loops} loops, {bhk_err_pct:.8f}% error)')
    print(f'  FD25 approximation:  sigma = 2**{log2(fd_sigma):.19f} ({fd_loops} loops, {fd_err_pct:.8f}% error)')
