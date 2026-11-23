# Reviewing Pull Requests

Anyone is welcome to review a *Landlab* pull request, not just core
maintainers. This page is a guide for what to look for and how to give
feedback that's useful to the author.

## Before you dig into the diff

* **Draft or ready?** PRs should start in draft. If it's still a draft,
  comment on direction.
* **Does it do one thing?** If a PR mixes unrelated changes, ask for it to
  be split before reviewing either half in detail.
* **Is the problem clear?** You should know *why* the change is needed, not
  just what it does. Ask if you don't.
* **Is CI green?** Don't spend time on a manual read of code that doesn't
  pass its own tests yet.
* **Is the history readable enough to review?** We use real merge
  commits, so branch commits become permanent history. Only ask for a
  rebase if messy history actually gets in the way of reviewing (e.g. a bug
  introduced then fixed in a later commit).

## Docstrings and tests

* Every public function/class/method has a *numpy* style docstring.
* Doctests are short and instructive, not just a correctness check.
* Unit tests cover every branch and error path, not just the happy path.
* Expected values in tests are worked out independently (e.g. by hand),
  not just pinned from the code's current output.
* Tests run on more than one grid type where relevant (`RasterModelGrid`
  and `HexModelGrid` at minimum).
* Function size is a judgment call but shorter is better.
  Flag it when a function bundles distinct concerns (e.g. validation, computation,
  I/O) so tightly that you can't test one without the others.
* New code comes with new tests, not just a coverage number that happens to
  hold.

## Component and utility API conventions

Check new or changed components against the component rules in the docs.

* `super().__init__(grid, ...)` is called.
* Public attributes are properties, not bare instance attributes.
* `_info` metadata is complete and dimensionally consistent with the rest
  of landlab.
* The component reads every input field and writes every output field it
  declares.
* Component methods follow an existing naming pattern (`run_one_step`,
  `update`, `calc_*`).

## Correctness

* **Conservation and symmetry.** Does a test actually check conserved
  quantities or expected symmetries, or just match a saved snapshot?
* **Boundary conditions and edge cases.** Closed boundaries, grids, all-zero
  input, etc.
* **Units and dimensional consistency**, especially across a field that will
  be passed between components.
* **Numerical stability**, does a scheme have an implicit stability
  assumption that should be checked or documented?
* **Reasonable defaults** for a first-time user, not just whatever the
  author last tested with.

Remember that it's fine for a review to say "I can't evaluate the numerics
here, can someone else look?" rather than approving on faith.

## Giving feedback

* Explain *why*, not just *what*.
* Flag a needed redesign early, before a line-by-line pass.
* Approve if you'd be fine seeing it merged as-is, use "request
  changes" for must-fix items.
