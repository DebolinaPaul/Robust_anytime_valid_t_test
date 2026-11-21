#!/usr/bin/env python3
"""Export shifted and one-sided reference values from the RPL38 Python engine.

Run from the eVille R-package root:
    python tools/export_shifted_reference.py --engine-dir /path/to/RPL38_final

The old engine's Prior is conditioned on theta>0. For a full prior with
location m, combine the positive-half mixture with the reflected mixture,
weighted by their UNCONDITIONED prior masses. For a theta<0 alternative,
reflect both t and m. Lebesgue halves are added without normalizing.
These are fixed synthetic cases, not a reanalysis of DepMap observations.
"""
from __future__ import annotations
import argparse
import csv
from functools import lru_cache
import hashlib
import importlib
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, log_ndtr
from scipy.stats import norm, cauchy


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def gaussian_radial_oracle(t, n, location, scale, alternative):
    """Independent positive radial integration after analytic Gaussian mixing.

    The effect-size integral is done analytically. The remaining chi_n
    integral is evaluated with adaptive quadrature, unlike engine.py's
    effect-size Gauss-Legendre rule. Truncation adds a normal CDF factor.
    """
    if alternative == 'less':
        return gaussian_radial_oracle(-t, n, -location, scale, 'greater')
    u = math.sqrt(n) * t / math.hypot(math.sqrt(n - 1), t)
    A = 1 + n * scale**2
    a = 1 - scale**2 * u**2 / A
    linear = location * u / (A * math.sqrt(a))
    constant = (-0.5 * math.log(A) - n * location**2 / (2*A)
                - n/2 * math.log(a) - (n/2-1)*math.log(2) - gammaln(n/2))
    if alternative == 'greater':
        constant -= float(log_ndtr(location / scale))
    def log_integrand(v):
        if v <= 0:
            return -np.inf
        out = constant + (n-1)*math.log(v) - v*v/2 + linear*v
        if alternative == 'greater':
            arg = (location + scale**2*u*v/math.sqrt(a)) / (scale*math.sqrt(A))
            out += float(log_ndtr(arg))
        return out
    raw_mode = (linear + math.hypot(linear, 2*math.sqrt(n-1)))/2
    bound = max(10., 3*raw_mode+10.)
    opt = minimize_scalar(lambda v: -log_integrand(v), bounds=(1e-12,bound),
                          method='bounded', options={'xatol':1e-12})
    mode = opt.x
    shift = log_integrand(mode)
    f = lambda v: math.exp(log_integrand(v)-shift)
    left, el = quad(f, 0., mode, epsabs=1e-12, epsrel=2e-12, limit=300)
    right, er = quad(f, mode, np.inf, epsabs=1e-12, epsrel=2e-12, limit=300)
    value = left+right
    if el+er > 1e-9 * value:
        raise ArithmeticError('Independent radial quadrature did not resolve the integral')
    return shift+math.log(value)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path,
                        default=Path('tests/testthat/fixtures'))
    args=parser.parse_args()
    root=args.engine_dir.expanduser().resolve(strict=True)
    source={name:root/name for name in ('engine.py','improper_flat.py')}
    for p in source.values():
        if not p.is_file():
            parser.error(f'Missing {p}')
    sys.path.insert(0,str(root))
    engine=importlib.import_module('engine')
    flat=importlib.import_module('improper_flat')
    if Path(engine.__file__).resolve()!=source['engine.py']:
        raise RuntimeError('Unexpected imported engine.py')
    if Path(flat.__file__).resolve()!=source['improper_flat.py']:
        raise RuntimeError('Unexpected imported improper_flat.py')
    out=args.output_dir.expanduser()
    out.mkdir(parents=True,exist_ok=True)

    @lru_cache(maxsize=None)
    def positive(t,n,family,m,s,order):
        prior=engine.Prior(family,m,s)
        return float(engine.log_mixture(t,n,prior,quadrature_order=order))

    def mixture(t,n,family,m,s,alternative,order):
        if family=='Lebesgue':
            plus=float(flat.log_flat_mixture(t,n))
            minus=float(flat.log_flat_mixture(-t,n))
            return {'greater':plus,'less':minus,
                    'two.sided':float(np.logaddexp(plus,minus))}[alternative]
        if alternative=='greater':
            return positive(t,n,family,m,s,order)
        if alternative=='less':
            return positive(-t,n,family,-m,s,order)
        dist=norm if family=='Gaussian' else cauchy
        log_pos=float(dist.logsf(-m/s))
        log_neg=float(dist.logcdf(-m/s))
        return float(np.logaddexp(log_pos+positive(t,n,family,m,s,order),
                                 log_neg+positive(-t,n,family,-m,s,order)))

    cases=[(n,t,family,m,s,alternative)
           for n in (2,5,12) for t in (-2.,0.,2.)
           for family in ('Gaussian','Cauchy')
           for m in (-3.,0.,3.) for s in (.7,1.4)
           for alternative in ('two.sided','greater','less')]
    cases += [(5,t,family,m,1.,alternative)
              for t in (-2.,2.) for family in ('Gaussian','Cauchy')
              for m in (-10.,10.)
              for alternative in ('two.sided','greater','less')]
    cases += [(n,t,'Lebesgue',0.,1.,alternative)
              for n in (2,5,12) for t in (-2.,0.,2.)
              for alternative in ('two.sided','greater','less')]
    rows,refinements,oracle_errors=[],[],[]
    for n,t,family,m,s,alternative in cases:
        low=mixture(t,n,family,m,s,alternative,32)
        high=mixture(t,n,family,m,s,alternative,64)
        refinements.append(abs(high-low))
        if family=='Gaussian':
            oracle=gaussian_radial_oracle(t,n,m,s,alternative)
            oracle_errors.append(abs(high-oracle))
        rows.append(dict(n=n,t=t,prior=family.lower(),location=m,scale=s,
                         alternative=alternative,log_m=high))
    max_refinement=max(refinements)
    max_oracle=max(oracle_errors)
    if max(max_refinement,max_oracle)>2e-8:
        raise ArithmeticError(f'Evidence check failed: refine={max_refinement}, oracle={max_oracle}')
    write_csv(out/'shifted_evidence.csv',rows)

    @lru_cache(maxsize=None)
    def calibrate(family,m,s,lc,order):
        return float(engine.calibration_direct(lc,engine.Prior(family,m,s),
                                               n0=2,quadrature_order=order))
    calibration_rows,cal_errors=[],[]
    for family in ('Gaussian','Cauchy'):
        for m in (-3.,0.,3.):
            for alternative in ('greater','less'):
                signed_m=m if alternative=='greater' else -m
                for lc in (-2.,0.,2.,math.log(20)):
                    low=calibrate(family,signed_m,1.4,lc,32)
                    high=calibrate(family,signed_m,1.4,lc,64)
                    cal_errors.append(abs(high-low))
                    calibration_rows.append(dict(prior=family.lower(),location=m,
                        scale=1.4,alternative=alternative,log_c=lc,g=high))
    for alternative in ('two.sided','greater','less'):
        for lc in (-2.,0.,2.,math.log(20)):
            if alternative=='two.sided':
                q=min(1.,math.sqrt(math.pi)*math.exp(-lc))
                value=(2/math.pi)*(math.asin(math.sqrt(q))+math.sqrt(q*(1-q)))
            else:
                value=float(flat.flat_calibration(lc))
            calibration_rows.append(dict(prior='lebesgue',location=0.,scale=1.,
                alternative=alternative,log_c=lc,g=value))
    if max(cal_errors)>2e-8:
        raise ArithmeticError(f'Calibration order agreement failed: {max(cal_errors)}')
    write_csv(out/'shifted_calibration.csv',calibration_rows)
    report={
        'scope':'shifted Gaussian/Cauchy; two-sided and conditioned one-sided mixtures; full/half Lebesgue',
        'source_files_sha256':{name:hashlib.sha256(path.read_bytes()).hexdigest()
                               for name,path in source.items()},
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'quadrature_orders':[32,64],
        'evidence_rows':len(rows),
        'calibration_rows':len(calibration_rows),
        'max_log_m_order_difference':max_refinement,
        'gaussian_independent_radial_checks':len(oracle_errors),
        'max_gaussian_radial_log_difference':max_oracle,
        'max_calibration_order_difference':max(cal_errors),
        'full_proper_mixture':'positive_prior_mass*Mplus(t;m) + negative_prior_mass*Mplus(-t;-m)',
        'less_mixture':'Mplus(-t;-m)',
        'lebesgue_two_sided':'Mplus(t)+Mplus(-t)',
        'calibration_note':'Proper calibration references only cover one-sided priors. They use the actual Python calibration_direct, a numerical approximation, not a certified bound. Full and half Lebesgue are analytic.'}
    (out/'shifted_reference_metadata.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
