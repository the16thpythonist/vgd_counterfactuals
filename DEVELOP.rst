Testing with Pytest
===================

During development it sometimes makes sense to run the tests only for a single version of python first. This can 
be done by using [pytest](https://docs.pytest.org/en/stable/) directly. To do this first make sure that pytest 
is installed in the current environmnent:

.. code-block:: bash

    uv tool install pytest

Then you can invoke the ``pytest`` command on the tests folder:

.. code-block:: bash

    pytest tests


Linting with Ruff
=================

Linters apply a set of rules on the codebase to check if the written code complies with certain coding guidelines. 
In general linters are used to keep and enforce a set of language-specific best practices across a code base and 
across a set of different developers.
This project uses the [ruff](https://github.com/astral-sh/ruff) linter, which can be installed like this:

.. code-block:: bash

    uv tool install ruff

To check the code against the linting rules use the ``ruff check`` command in the top-level folder:

.. code-block:: bash

    ruff check .


Bumping Version for a new Release
================================= 

To release a new version of the package, the version string has to be updated throughout all the different 
places where this version string is used in the text. In this project, this is handled automatically 
using the [bump-my-version](https://github.com/callowayproject/bump-my-version) tool, which can be 
installed like this:

.. code-block:: bash

    uv tool install bump-my-version

One of the following commands can then be used to bump the version either for a patch, minor or major release: 

.. code-block:: bash

    bump-my-version bump -v patch
    bump-my-version bump -v minor
    bump-my-version bump -v major

The configuration of which files are being updated and how the version is parsed etc. can be found in a 
tool section of the ``pyproject.toml``


Building a new Package Version
==============================

Before a new version of the package can be published on PyPi for example, the code has to be built first. This 
can be done with uv's ``build`` command like this:

.. code-block:: bash

    uv build --python=3.10

If it doesn't already exist, this command will create a new ``dist`` folder where the built tarball and wheel of 
the current version (as defined in the pyproject.toml file) are saved.


Publishing a new Version to PyPi
================================

[twine](https://twine.readthedocs.io/en/stable/) is a python library that is specifically intended for publishing python 
packages to the package indices such as PyPi. Twine can be installed like this:

.. code-block:: bash

    uv tool install twine

After this the ``twine`` command is available:

.. code-block:: bash

    twine --help

**Checking the distribution. ** Twine assumes that the built distribution files (tarball and wheel) already exist in the 
project's ``dist`` folder (see "Building a New Package Version"). The ``twine check`` command can be used to check 
these distribution files for correctness before actually uploading them. This command will for example check the 
syntax of the README file to make sure it can be properly rendered on the PyPi website.

.. code-block:: bash

    twine check dist/*
    
**Uploading to PyPi.** Finally, the ``twine upload`` command can be used to actually upload the distribution files 
to the package index.

    twine upload --username='__token__' --password='[your password]' dist/*


Documentation with MkDocs
=========================

The documentation is done with [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/). The documentation configuration 
can be found in the ``mkdocs.yml`` file and the actual markdown files are in the top-level ``docs`` folder of the project.

**Local Development.** To view the local dev version of the documentation, you can use the ``mkdocs serve`` command:
like this:

.. code-block:: bash

    mkdocs serve    

This will start a development web server to serve the static doc files which can then be viewed with a browser.

**Publishing to Github Pages.** The production version of the documentation is hosted on Github Pages. Once a sufficient update 
of the documentation was written locally, these changes can be published to the Gh Pages branch of the remote repository by 
using the following command:

.. code-block:: bash

    mkdocs gh-deploy --force