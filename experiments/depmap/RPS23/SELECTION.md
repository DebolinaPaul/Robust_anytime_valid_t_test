# How this example was selected

This is a retrospective search for a second illustration, as requested. It
must not be described as a prespecified test of RPS23 selected independently
of the target data.

1. Start with the complete DepMap Public 24Q4 gene-effect matrix and the same
   Model.csv metadata as the supplied RPL38 package. Consider recorded Female
   and Male cohorts separately; keep one ModelID per PatientID. Retain lineages
   with at least ten observations. For this preliminary screen, restrict to
   genes with complete finite scores on the relevant sex-specific models.
2. Find a weakest standardized target effect between -4.5 and -1.25, Gaussian
   source location below -4, all other source effects below -2, and a positive
   source-target location-scale discrepancy above 2.25. Exclude RPL38. This
   produced 243 gene/sex/target candidates.
3. Evaluate 43 initial shortlisted combinations with the unchanged Gaussian
   evidence/calibration functions, Rouder and Lebesgue benchmarks, and the
   same 100 paired random orders. Then evaluate 25 additional combinations
   with Gaussian source location below -6.5. The seeds and ascending ModelID
   illustration order were not tuned.
4. Check all five fitted priors for promising candidates. Select RPS23/Male,
   with Prostate as the illustrative target, because it has substantial
   standardized-effect conflict, earlier mean eVille stopping for every
   fitted prior, an interior horseshoe scale for that target, and separated
   Gaussian crossings in the existing fixed order.
5. Run the full, unchanged leave-one-lineage-out analysis for this gene:
   all 18 eligible lineages, all seven priors/benchmarks, and 100 paired random
   orders per lineage. Do not drop eligible lineages based on stopping times.

The source priors in each final fold exclude that target lineage. This avoids
target leakage into the fitted prior for a fixed fold, but does not undo the
retrospective selection of gene and illustrative lineage. Nominal per-test
calibration is not a claim of selective or family-wise error control for the
search. The included figures describe the selected observed data, not a
controlled simulation of increasing conflict.
