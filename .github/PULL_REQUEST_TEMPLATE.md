# Pull Request Guidelines

Thank you for contributing! Please follow these guidelines to help keep our
codebase clean, consistent, and maintainable.

## What Makes a Good Pull Request

- **Small and focused**: Address one issue or feature per pull request. Avoid
  bundling multiple unrelated changes.
- **Clear purpose**: The pull request should explain *why* a change is needed,
  not just *what* was changed.
- **Draft first**: Open your pull request in **draft mode** so others can
  follow progress and provide early feedback.
- **Readable history**: Prefer a clean commit history. Consider squashing
  or rebasing before final review.
- **Well-tested and documented**: All code changes should be accompanied
  by appropriate tests and updated documentation.

---

## Pull Request Checklist

Please confirm the following before marking the pull request as "Ready for Review"
(if any of the following items aren't relevant for your contribution please still
tick them so we know you've gone through the checklist):

- [ ] The PR is opened in **draft mode**
- [ ] The PR addresses a single feature, fix, or issue
- [ ] Code has been linted and is free of style issues (`nox -s lint`)
- [ ] All tests pass locally (`pytest`, or `nox -s test`)
- [ ] Documentation builds without errors (`nox -s docs-build`)
- [ ] A [**news fragment**](https://landlab.csdms.io/development/contribution/#news-entries)
      describing the change has been added
- [ ] Relevant docstrings, examples, or changelog entries have been updated
- [ ] Related issues or discussions are linked in the description (e.g., `Closes #1973`)

---

## Description

<!--
Provide a short summary of the change. Why is it needed? What does it do?
-->

---

## Related Issues

<!--
List any related issues, discussions, or pull requests.
Use keywords like "Closes #1973" to automatically close issues when merged.
-->
