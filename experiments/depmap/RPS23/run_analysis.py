#!/usr/bin/env python3
"""One final RPS23 analysis: source-only MLE fits and all held-out lineages.

Run: python -u run_analysis.py --permutations 100
The horseshoe uses source-only noncentral-t marginal maximum likelihood.
All outputs are consolidated; no per-lineage output directories are created.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from cohort import load_gene_data, partition, summarize, assess_cohort
from engine import log_mixture
from fitting import build_priors, IMPROPER_LABEL
from sequential import (
    lineage_seed, permutation_t_statistics, cached_log_threshold,
    boundaries_for_evaluated_orders, first_crossings,
    rejection_summary, paired_summary, DirectCalibration,
)


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def save_json(path, value):
    Path(path).write_text(json.dumps(
        json_safe(value), indent=2, allow_nan=False) + '\n')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def outcome_fractions(e, v, N):
    names = ['outcome_same_rejection_fraction', 'outcome_save1_fraction',
             'outcome_save2_fraction', 'outcome_save3plus_fraction']
    if v is None:
        return dict.fromkeys(names)
    both = (e <= N) & (v <= N)
    difference = v - e
    return dict(zip(names, [float(np.mean(both & (difference == k)))
                            for k in (0, 1, 2)]
                    + [float(np.mean(both & (difference >= 3)))]))


def analyse(args, out):
    data, audit = load_gene_data(
        args.gene_effect, args.models, args.gene, args.sex)
    available = sorted(data.OncotreeLineage.unique())
    cohort = assess_cohort(
        data, available, args.min_source_n, args.min_target_n)
    estimates = pd.DataFrame([
        {'target_lineage': lineage, **summarize(group.Chronos)}
        for lineage, group in data.groupby('OncotreeLineage', sort=True)
    ]).drop(columns='n')
    cohort = cohort.merge(estimates, on='target_lineage',
                          validate='one_to_one')
    cohort['cohort'] = np.where(
        cohort.target_eligible,
        np.where(cohort.target_n >= args.main_min_n, 'main', 'small'),
        'not_testable',
    )
    targets = cohort.loc[cohort.target_eligible, 'target_lineage'].tolist()
    if not targets:
        raise ValueError('No eligible target lineages')
    if args.example_target not in targets:
        raise ValueError('--example-target must be an eligible target lineage')

    settings = dict(
        gene=args.gene, sex=args.sex, alpha=args.alpha, n0=2,
        permutations=args.permutations, seed=args.seed,
        min_source_n=args.min_source_n, main_min_n=args.main_min_n,
        min_target_n=args.min_target_n, quadrature_order=args.quadrature_order,
        validation_orders=args.validation_orders, example_target=args.example_target,
        source_fitting=(
            "four full-family MLE fits; "
            "horseshoe constrained marginal MLE with noncentral-t sampling"
        ),
        target_order='ascending ModelID is fixed order 0; summaries use only random orders 1,...,B',
        permutation_seed_rule='SHA256(base seed, gene, sex, lineage); first 8 bytes little endian',
        observations_retained=len(data), tested_lineages=targets,
        main_lineages=cohort.loc[cohort.cohort.eq(
            'main'), 'target_lineage'].tolist(),
        small_lineages=cohort.loc[cohort.cohort.eq(
            'small'), 'target_lineage'].tolist(),
        input_sha256={'gene_effect': sha256(
            args.gene_effect), 'models': sha256(args.models)},
        conflict_definition='abs(target standardized sample mean - full-family prior location) / prior scale; descriptive, not a p-value; scale is not necessarily a standard deviation',
        interpretation='Permutation summaries condition on the observed sample; they are not repeated-sample power. True target effects are unknown. Main and small-lineage cohorts must be distinguished.',
        package='RPS23_final',
    )
    save_json(out / 'settings.json', settings)
    save_json(out / 'data_filter_audit.json', audit)
    data.to_csv(out / 'observations.csv.gz', index=False, compression='gzip')

    all_stops, summaries, fits, boundaries, example_paths = [], [], [], [], []
    numerical_checks = []
    for position, lineage in enumerate(targets, 1):
        print(
            f'[{position}/{len(targets)}] {lineage}: fitting source priors', flush=True)
        target, sources, source_summary, _ = partition(
            data, lineage, args.min_source_n)
        x = -target.Chronos.to_numpy(float)
        N = len(x)
        seed = lineage_seed(args.seed, args.gene, args.sex, lineage)
        _, ts, ns = permutation_t_statistics(x, args.permutations, seed)
        target_info = cohort.set_index('target_lineage').loc[lineage]
        priors, diagnostics = build_priors(
            -source_summary.theta_hat.to_numpy(float),
            source_summary.n.to_numpy(int),
        )
        diagnostic_by_prior = {record['prior']
            : record for record in diagnostics}
        local_summaries = []
        local_stops = {}

        for label, prior in priors.items():
            proper = label != IMPROPER_LABEL
            te = cached_log_threshold(
                prior, args.alpha, 2, args.quadrature_order)
            tv = -math.log(args.alpha) if proper else None
            if proper and te > tv + 1e-8:
                raise ArithmeticError(
                    'eVille threshold exceeds Ville threshold')
            be, bv = boundaries_for_evaluated_orders(
                ts, ns, prior, te, tv, args.quadrature_order)
            se = first_crossings(ts, be, N)
            sv = first_crossings(ts, bv, N) if proper else None
            if proper and np.any(se > sv):
                raise ArithmeticError('eVille crossing is later than Ville')
            local_stops[label] = (se, sv)

            # Check root-based decisions with direct likelihood-mixture values.
            for order in range(min(args.permutations, args.validation_orders) + 1):
                logs = np.array([log_mixture(float(t), int(n), prior,
                                             args.quadrature_order)
                                 for t, n in zip(ts[order], ns)])
                if not np.isfinite(logs).all():
                    raise ArithmeticError(
                        'Nonfinite evidence in direct boundary validation')
                for method, threshold, stops in [('eVille', te, se)] + (
                        [('Ville', tv, sv)] if proper else []):
                    hits = np.flatnonzero(logs >= threshold)
                    direct = int(ns[hits[0]]) if len(hits) else N + 1
                    if direct != int(stops[order]):
                        raise ArithmeticError(
                            f'Boundary mismatch: {lineage}/{label}/{order}/{method}')
                    numerical_checks.append(dict(target_lineage=lineage, prior=label,
                                                 permutation=order, method=method,
                                                 direct_and_boundary_agree=True))
                if order == 0 and lineage == args.example_target:
                    maxima = np.maximum.accumulate(logs)
                    pe = DirectCalibration(
                        prior, quadrature_order=args.quadrature_order).pvalue(maxima)
                    pv = np.exp(np.minimum(0., -maxima)
                                ) if proper else np.full_like(maxima, np.nan)
                    for j, n in enumerate(ns):
                        example_paths.append(dict(prior=label, n=int(n),
                                                  t_raw_chronos=-float(ts[0, j]), log_evidence=logs[j],
                                                  log_running_evidence=maxima[j], p_eville=pe[j], p_ville=pv[j]))

            for method, values in [('eVille', be)] + ([('Ville', bv)] if proper else []):
                for n, value in zip(ns, values):
                    status = ('not_needed_after_all_orders_rejected' if np.isnan(value)
                              else 'finite_root' if np.isfinite(value)
                              else 'above_searched_range' if value > 0 else 'below_searched_range')
                    boundaries.append(dict(target_lineage=lineage, prior=label, method=method,
                                           n=int(n), critical_t=value, boundary_status=status))

            base = dict(gene=args.gene, sex=args.sex, target_lineage=lineage, target_n=N,
                        cohort=target_info.cohort, source_n=len(sources),
                        source_lineages=len(source_summary), prior=label)
            fit_record = {**base, **diagnostic_by_prior[label], 'family': prior.family,
                          'is_proper': proper, 'full_prior_location_raw_theta': -prior.location if proper else None,
                          'scale': prior.scale if proper else None, 'log_eville_threshold': te,
                          'eville_threshold': math.exp(te), 'ville_threshold': 1 / args.alpha if proper else None}
            fits.append(fit_record)
            for order, e in enumerate(se):
                v = sv[order] if proper else None
                all_stops.append({**base, 'permutation': order,
                                  'eville_rejection': int(e) if e <= N else None,
                                  'ville_rejection': int(v) if proper and v <= N else None,
                                  'eville_rejected': bool(e <= N),
                                  'ville_rejected': bool(v <= N) if proper else None})

            e_summary = rejection_summary(se[1:], N)
            v_summary = rejection_summary(
                sv[1:], N) if proper else dict.fromkeys(e_summary)
            signed = (target_info.theta_hat + prior.location) / \
                prior.scale if proper else None
            local_summaries.append({**base,
                                    'theta_hat': target_info.theta_hat,
                                    'theta_ci95_low': target_info.theta_ci95_low,
                                    'theta_ci95_high': target_info.theta_ci95_high,
                                    **{'eville_' + k: v for k, v in e_summary.items()},
                                    **{'ville_' + k: v for k, v in v_summary.items()},
                                    **paired_summary(se[1:], sv[1:] if proper else None, N),
                                    **outcome_fractions(se[1:], sv[1:] if proper else None, N),
                                    'conflict_signed_location_scale': signed,
                                    'conflict_abs_location_scale': abs(signed) if proper else None,
                                    'original_eville_rejection': int(se[0]) if se[0] <= N else None,
                                    'original_ville_rejection': int(sv[0]) if proper and sv[0] <= N else None,
                                    'fit_method': diagnostic_by_prior[label]['actual_method'],
                                    })

        flat = local_stops[IMPROPER_LABEL][0][1:]
        flat_mean = float(np.minimum(flat, N).mean())
        for row in local_summaries:
            row['flat_mean_capped_at_N'] = flat_mean
            row['eville_mean_minus_flat'] = row['eville_mean_capped_at_N'] - flat_mean
            row['ville_mean_minus_flat'] = (row['ville_mean_capped_at_N'] - flat_mean
                                            if row['ville_mean_capped_at_N'] is not None else None)
            summaries.append(row)

    cohort.to_csv(out / 'cohort.csv', index=False)
    pd.DataFrame(summaries).to_csv(out / 'summary.csv', index=False)
    pd.DataFrame(fits).to_csv(out / 'prior_fits.csv', index=False)
    pd.DataFrame(all_stops).to_csv(
        out / 'stopping_times.csv.gz', index=False, compression='gzip')
    pd.DataFrame(boundaries).to_csv(
        out / 'critical_boundaries.csv.gz', index=False, compression='gzip')
    pd.DataFrame(example_paths).to_csv(out / 'example_paths.csv', index=False)
    save_json(out / 'boundary_checks.json', dict(checks=len(numerical_checks),
                                                 all_agree=True, checked_orders='original order plus first validation_orders random orders',
                                                 records=numerical_checks))
    return settings


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gene', default='RPS23')
    parser.add_argument('--sex', choices=['Female', 'Male'], default='Male')
    parser.add_argument('--gene-effect', type=Path,
                        default=root / 'data/CRISPRGeneEffect.csv')
    parser.add_argument('--models', type=Path, default=root / 'data/Model.csv')
    parser.add_argument('--out', type=Path, default=root / 'results')
    parser.add_argument('--permutations', type=int, default=100)
    parser.add_argument('--seed', type=int, default=20261009)
    parser.add_argument('--min-source-n', type=int, default=10)
    parser.add_argument('--main-min-n', type=int, default=10)
    parser.add_argument('--min-target-n', type=int, default=10)
    parser.add_argument('--example-target', default='Prostate')
    parser.add_argument('--alpha', type=float, default=.001)
    parser.add_argument('--quadrature-order', type=int, default=48)
    parser.add_argument('--validation-orders', type=int, default=2)
    parser.add_argument('--no-figures', action='store_true')
    parser.add_argument('--overwrite', action='store_true',
                        help='Replace only a recognized completed RPS23_final output folder.')
    args = parser.parse_args()
    if (args.permutations < 1 or not 0 < args.alpha < 1 or args.seed < 0
            or min(args.min_source_n, args.min_target_n, args.main_min_n) < 10
            or args.validation_orders < 0 or args.quadrature_order < 24):
        parser.error(
            'Require positive permutations, alpha in (0,1), nonnegative seed/validation-orders, sizes >=10 and quadrature-order >=24')
    out = args.out.resolve()
    if out.exists() and any(out.iterdir()):
        if not args.overwrite:
            parser.error(
                'Output already contains files; choose a new --out or use --overwrite for this package\'s completed output')
        metadata = out / 'settings.json'
        manifest = out / 'manifest.json'
        if (not metadata.exists() or not manifest.exists()
                or json.loads(metadata.read_text()).get('package') != 'RPS23_final'
                or json.loads(manifest.read_text()).get('complete') is not True):
            parser.error('Refusing to overwrite an unrecognized directory')
    if out == root or root.is_relative_to(out):
        parser.error('Output must not be the code directory or its parent')
    out.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='rps23_run_', dir=out.parent))
    try:
        analyse(args, staging)
        from validate import validate_results
        validation = validate_results(staging)
        save_json(staging / 'validation.json', validation)
        if not args.no_figures:
            from figures import make_figures
            make_figures(staging)
        save_json(staging / 'manifest.json', {
            'complete': True,
            'files': sorted(str(p.relative_to(staging)) for p in staging.rglob('*') if p.is_file()),
        })
        if out.exists():
            shutil.rmtree(out)
        staging.replace(out)
    except Exception:
        print(
            f'Incomplete run retained for diagnosis: {staging}', file=sys.stderr)
        raise
    print(f'Completed and verified: {out}', flush=True)


if __name__ == '__main__':
    main()
