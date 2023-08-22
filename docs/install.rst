Installation
============

From the Python Package Index (PyPI)
------------------------------------

Probably the easiest way to get up and running::

    pip install fair


From anaconda
-------------

    conda install -c chrisroadmap fair


From GitHub
-----------

Download and install
~~~~~~~~~~~~~~~~~~~~

The latest release can be obtained from https://github.com/OMS-NetZero/FAIR/releases as zip or tarball files, or the most current unreleased version can be cloned from https://github.com/OMS-NetZero/FAIR.

Developing
~~~~~~~~~~

1. Fork the repository, then clone it to your local disk. `cd` to the `FAIR` directory in your working copy.
2. Create a new branch for your changes::

    git checkout -b <branch_name>

3. Optional, but we highly recommend developing FaIR in a virtual environment to keep your base installation of `python` nice and clean.
4. Install `fair` in development mode::

    pip install -e .[dev]

5. Make your changes.
6. Write a test that tests your new feature (in the ``tests`` directory of the repository).
7. Format your code, and run tests locally::

    make format
    make checks
    make tests
    make test_notebooks

If you find errors at this point, they will need fixing before GitHub will allow merging to the `master` branch. Running the test suite ensures that your code change does not break or change existing functionality.

8. Commit and push your changes::

    git add <file>
    git commit -m "informative commit message"
    git push origin <branch_name>

9. Create a pull request on the parent repository to merge in your changes. You will see there's a checklist of processes to go through...
10. One of which is adding a line to `CHANGELOG.rst` with a summary of the changes made.
11. The checks and tests will run using GitHub actions. If all pass, you should be able to submit the pull request for review.
12. If the codeowners are happy, the branch will be merged in.

TODO: Check out the (currently non-existent) contributing guide, but it's basically this.
