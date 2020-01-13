.. _joss_workflow:

===============================================================================
Workflow to submit your Landlab contribution to Journal of Open Source Software
===============================================================================

The `Journal of Open Source Software (JOSS) <https://joss.theoj.org>`_ is a
venue for publications on open source software. If you are considering
preparing a contribution to Landlab for a JOSS publication, please make sure to
review the JOSS `submission requirements
<https://joss.readthedocs.io/en/latest/submitting.html#submission-requirements>`_


Parts of Landlab have already been published by JOSS, and multiple parts of
Landlab may be under consideration at JOSS at the same time. Under these
circumstances there may be conflicts with the required `paper.md` file.

For this reason, the workflow below has been developed for those who choose to
submit their Landlab contribution (e.g., a Landlab component, utility) to this
journal.


- Create a pull request to merge the contribution into Landlab. See
  :ref:`Develop your own component or utility <landlab_component_dev_page>`.
- Following the above merge, create a new branch for the JOSS review process.
  Name the paper ``paper.md``, which follows JOSS requirements, and place this
  file in a new folder, ``/landlab/joss/<name_of_contribution>/`` with your other
  JOSS submission materials (e.g., bib, figures if any).
- Submit to JOSS. The submitting author must comment on the pre-review issue
  that the paper is on a branch, and they should build the branch with this
  command

  .. code-block:: bash

      $ @whedon prepare pdf from branch my-custom-branch-name

- Once the JOSS review is complete and if the paper is accepted:
- Create a pull request to merge the JOSS review branch into the ``master``
  branch.
- Merge ``master`` into ``release``, and update version number. See
  :ref:`Create a Landlab release <dev_releases>`.
- Ensure that the release is correctly created and distributed.
- Update the JOSS review pull request with the correct archive DOI and
  version number.
- JOSS publishes the paper.
- Once publication is complete, make another pull request to move and rename
  the ``paper.md`` to ``landlab/joss/published/<name_of_module_or_component>/<new_paper_name>.md``.
  The paper can no longer be named ``paper.md`` because JOSS expects paper in
  review to have this name, and to avoid conflicts with other Landlab paper
  contributions in review. Note that the paper folder should be moved to
  within the ``published`` directory.
