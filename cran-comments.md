## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE

      Maintainer: 'Olatunji Johnson <teejaystatistics@gmail.com>'
      New submission

  This is the expected "New submission" note for a first CRAN release; there are no
  other notes. The package URL and BugReports
  (https://github.com/olatunjijohnson/spLeakage) and the documentation website
  (https://olatunjijohnson.github.io/spLeakage/) are live and resolve.

## Test environment

* local: R 4.5.2 on Linux

## Notes

* The optimism emulator shipped in `R/sysdata.rda` is a compact (~85 KB) gradient-
  boosted model serialised as raw bytes; `xgboost` is used only at predict time and
  is declared in Suggests (the function errors informatively if it is unavailable).
* Heavier dependencies used only to build bundled artefacts and case studies
  (`ranger`, `mgcv`, `fields`, `xgboost`, `malariaAtlas`) are in Suggests or used
  only under `data-raw/`.
