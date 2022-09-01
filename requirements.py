#! /usr/bin/env python
import argparse
import os


def requirements(extras, include_required=False):

    sections = Requirements().get_sections(
        extras, include_required=not extras or include_required
    )

    print("# Requirements extracted from pyproject.toml")
    for section, packages in sections.items():
        print(f"# {section}")
        print(os.linesep.join(sorted(packages)))


class Requirements:
    def __init__(self):
        tomllib = _find_tomllib()

        with open("pyproject.toml", "rb") as fp:
            project = tomllib.load(fp)["project"]
        self._dependencies = Requirements._clean_dependency_names(
            project.get("dependencies", [])
        )
        self._optional_dependencies = {
            name: Requirements._clean_dependency_names(packages)
            for name, packages in project.get("optional-dependencies", {}).items()
        }

    @staticmethod
    def _clean_dependency_names(packages):
        return [package.split(";")[0].strip() for package in packages]

    @property
    def required(self):
        return tuple(self._dependencies)

    def get_optional(self, extra):
        return tuple(self._optional_dependencies[extra])

    def get_sections(self, sections=None, include_required=False):
        sections = sections or []
        reqs = {}
        if include_required:
            reqs["[project.dependencies]"] = self.required
        for section in sections:
            reqs[f"[project.optional-dependencies] {section}"] = self.get_optional(
                section
            )
        return reqs


def _find_tomllib():
    try:
        import tomllib
    except ModuleNotFoundError:
        import tomli as tomllib
    return tomllib


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract requirements information from pyproject.toml"
    )
    parser.add_argument("extras", type=str, nargs="*")
    args = parser.parse_args()

    requirements(args.extras)
