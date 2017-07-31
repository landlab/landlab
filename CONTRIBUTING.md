# Contributing to landlab

Thank you for contributing to the landlab! We appreciate
your help as this is largely as volunteer effort! :heart: :heart: :heart:

# How to contribute

## Reporting Bugs

Before creating a bug report, please do at least a cursory check that the
bug has not already been reported. If it has, add a comment to the existing
issue instead of opening a new one.

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
*  **Provide a [minimal example](https://stackoverflow.com/help/mcve)
   that demonstrates the steps** as, for example, a bash script
   along with input files. This example should reproduce your
   problem with as few lines of code as possible and easily
   reproducible my another person.
*  **Describe the behavior you are seeing after these steps**.
*  **Describe the behavior you expect to see after these steps**.

Additionally, the answers to the following questions about your run
environment will be helpful.

*  **Which version of landlab are you using?** This could be a specific
   git sha or a release number.
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
   and solution.  Include the relevant issue number.

## Styleguides

Use the [PEP8 style guide](https://www.python.org/dev/peps/pep-0008/).
You may want to use tools like
[Prospector](http://prospector.landscape.io/en/master/) to help out
with this.

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

## Adding new components

If you would like to create a new component, we have just a few
conventions that we would like you to follow.

* Create a new folder under `landlab/components` that will hold your
  new component.
* Be aware that your component will have to be able to work with both
  Python 2 and Python 3.

### Python 2 vs. Python 3

* We use `six` to maintain a single code base that will work with
  both Python 2 and 3.
* All tests will be run with Python 2.7, 3.5, and 3.6.
* If you need to introduce a new dependency, that dependency must
  be compatible with Python 2 and 3.

### Documentation

* All public functions, classes, methods, etc. must have a docstring
  that follows the [numpydoc](https://github.com/numpydoc/numpydoc)
  conventions.
* Every `.py` file must contain a module-level docstring that the top
  of the file that describes what the purpose of the file is.

### Testing

* All contributed code must be well tested. This should be done through
  both doctests as well as more standard unit tests through `nose`.
* Doctests should be short, easy-to-read tests that are instructive
  to a user.
* Unit tests should be significanly more extensive and give your
  new code thorough testing.


Thanks! :heart: :heart: :heart:

The landlab team
