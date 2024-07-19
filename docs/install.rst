Installation
============

`fair` requires `python` 3.8-3.12.

From anaconda (recommended)
---------------------------

::

    conda install -c conda-forge fair


From the Python Package Index (PyPI)
------------------------------------

::

    pip install fair


From GitHub
-----------

The latest release can be obtained from https://github.com/OMS-NetZero/FAIR/releases as zip or tarball files, or the most current unreleased version can be cloned from https://github.com/OMS-NetZero/FAIR.

To install::

    pip install -e .[dev]

This will also give access to the notebooks in the `examples` directory, so you can experiment with the model using these as a starting point. For this to work, you should launch notebooks from the same directory in which you install `FAIR` locally.

Guide for developers
--------------------

1. Navigate to the GitHub repository, fork to your personal account, then clone it to your local disk. `cd` to the `FAIR` directory in your working copy.
2. Create a new branch for your changes::

    git checkout -b <branch_name>

3. The development package includes a Makefile which will set up a virtual environment for you along with other tools, keeping your base installation clean. `python` 3.11 is used for developing FaIR. To set this up, run::

    make venv

4. Make your code changes.
5. Write a test that tests your new feature (in the ``tests`` directory of the repository).
6. Format your code, and run tests locally::

    make format
    make checks
    make test
    make test_notebooks

If you find errors at this point, they will need fixing before GitHub will allow merging to the `master` branch. Running the test suite ensures that your code change does not break or change existing functionality. The remote and local tests will also fail if you have not increased the code coverage (basically, if you haven't written a test for your change).

7. Commit and push your changes::

    git add <file>
    git commit -m "informative commit message"
    git push origin <branch_name>

8. Create a pull request on the parent repository to merge in your changes. You will see there's a checklist of processes to go through...
9. One of which is adding a line to `CHANGELOG.rst` with a summary of the changes made.
10. The checks and tests will run using GitHub actions. If all pass, you should be able to submit the pull request for review.
11. If the codeowners are happy, the branch will be merged in.
