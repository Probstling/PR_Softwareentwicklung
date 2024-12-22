.. PropScript documentation master file, created by
   sphinx-quickstart on Sat Dec 21 18:16:50 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PropScript documentation
========================

Add your content using ``reStructuredText`` syntax. See the
`reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
documentation for details.


Overview
--------

This project contains three main Python scripts:

- **extract.py**: Extracts sequence specific properties and lists them in a pandas DataFrame.
- **figures.py**: Creates violin plots from the extracted data.
- **features.py**: Maps extracted features to the DataFrame and exports it as a CSV file.

Scripts
-------

.. toctree::
   :maxdepth: 1
   :caption: Scripts:

   extract
   figures
   features

Script Details
--------------

.. literalinclude:: ../extract.py
   :language: python
   :caption: "extract.py"

.. literalinclude:: ../figures.py
   :language: python
   :caption: "figures.py"

.. literalinclude:: ../features.py
   :language: python
   :caption: "features.py"