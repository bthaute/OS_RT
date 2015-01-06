.. pycontroltools__sphinx_README-label:

*****************************************************
How to build the Sphinx documentation for pycontroltools
*****************************************************

The documentation for the **pycontroltools** module is created using Sphinx.
To be able to use the tools Sphinx has to be installed
(http://sphinx-doc.org/latest/install.html)
and initialized.

Initializing Sphinx
---------------------
To initialize the Sphinx folder hierarchy for the project,
you preferably create a new folder for the documentation in the main folder of
your project.
(The script files containing the docstrings have to be in the same directory
as the documentation folder or in the documentation folder itself.)

Now start a command console in that folder and type:

.. code-block:: bash

    sphinx-quickstart
	
Build HTML- and LaTeX-Files
--------------------------------

Copy and replace the files from the **pycontroltools_sphinx** folder of the
repositorie into your documentation folder.

To build the HTML-files type:

.. code-block:: bash

    make html
	
The build-files appear in the folder **_build/html**.

To build the LaTeX-Files type:

.. code-block:: bash
    
    make latex
	
The build-files appear in the folder **_build/latex**.
If a LaTeX distribution is installed on your system, you can now use
the LaTeX files to build (and edit) the document with your own LaTeX editor.