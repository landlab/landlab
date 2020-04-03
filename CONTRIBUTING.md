# Contributing to landlab

Thank you for contributing to Landlab! We appreciate
your help as this is largely as volunteer effort! :heart: :heart: :heart:

# How to contribute

## Reporting Bugs

Before creating a bug report, please do at least a cursory check that the
bug has not already been reported by searching the Issues portion of the
GitHub repository. If it has, add a comment to the existing issue instead of
opening a new one.

### Submitting a Bug Report

Bugs are tracked as
[GitHub issues](https://guides.github.com/features/issues/). After you've
determined you've found a new bug, please open a
[new issue](https://github.com/landlab/landlab/issues).

Explain the problem and include additional details to help maintainers
reproduce the problem. Here are some items that will make it easier
to track down the source of the problem.

*  **Use a clear and descriptive title** for the issue that identifies the
   problem.
*  **Describe the exact steps that reproduce the problem**.
*  **Provide a [minimal example](https://stackoverflow.com/help/minimal-reproducible-example)
   that demonstrates the steps** as, for example, a bash script
   along with input files. This example should reproduce your
   problem with as few lines of code as possible and easily
   reproducible my another person. Such an example almost certainly will not
   include an input file or any dependencies beyond those required by the
   `landlab_dev` conda environment.
*  **Describe the behavior you are seeing after these steps**.
*  **Describe the behavior you expect to see after these steps**.

Additionally, the answers to the following questions about your run
environment will be helpful.

*  **Which version of landlab are you using?** This could be a specific
   git sha or a release number. The best way to find this information is to
   import landlab and evaluate `landlab.__version__`
*  **What is he name and version of you OS?**
*  **What compiler are you using?**
*  **How did you build landlab (if using the development version)?**


## Submitting Changes

:tada: Whoa! This is great! We love it when folks contibute code! :tada:

Changes to landlab should be submitted as
[pull requests](http://help.github.com/pull-requests/)).

*  Create a GitHub issue that describes what you propose to do.
*  Create a topic branch that contains your changes.
*  Open a new [GitHub pull request](https://github.com/landlab/landlab/pull/new/master).
*  Ensure the pull request description clearly describes the problem
   and solution. Include the relevant issue number.

### Git Commit Messages

* Use the present tense ("Add feature" not "Added feature")
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally
* For fun, consider starting the commit message with an applicable emoji:
    * :art: `:art:` when improving the format/structure of the code
    * :racehorse: `:racehorse:` when improving performance
    * :non-potable_water: `:non-potable_water:` when plugging memory leaks
    * :memo: `:memo:` when writing docs
    * :penguin: `:penguin:` when fixing something on Linux
    * :apple: `:apple:` when fixing something on macOS
    * :checkered_flag: `:checkered_flag:` when fixing something on Windows
    * :bug: `:bug:` when fixing a bug
    * :fire: `:fire:` when removing code or files
    * :green_heart: `:green_heart:` when fixing the CI build
    * :white_check_mark: `:white_check_mark:` when adding tests
    * :shirt: `:shirt:` when removing linter warnings

### Pull Request Messages

  * Rename the pull request and provide a comment that synthesizes what
    the pull request changes or adds. This helps us synthesize what
    changes have occured between Landlab releases.

## Adding new components

If you would like to create a new component, we a few conventions that we would
like you to follow.

Please visit [this part](https://landlab.readthedocs.io/en/master/development/index.html)
of the main Landlab documentation page to read about developer installation,
guidelines to contributing code, and our software development practices.

**Landlab 2 is Python >=3.6 only.**

Thanks! :heart: :heart: :heart:

The Landlab team
