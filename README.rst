==================================
Schrodinger Maestro scipion plugin
==================================

Programs coming from Schrodinger

===================
Install this plugin
===================

You will need to use `3.0.0 <https://github.com/I2PC/scipion/releases/tag/v3.0>`_ version of Scipion to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block:: 

      scipion installp -p scipion-chem-schrodinger
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. Download repository: 

.. code-block::

            git clone https://github.com/scipion-chem/scipion-chem-schrodinger.git

2. Install:

.. code-block::

            scipion installp -p path_to_scipion-chem-schrodinger --devel

- **Binary files** 

Atom_struct_utils plugin is a pure Python module, no binary files are required.

The **Schrodinger software must be installed separately**, as it requires a license.
You must either install it in the EM_ROOT (typically: SCIPION_HOME/software/em/Schrodinger2021-3)
or define the Schrodinger path in scipion.conf as SCHRODINGER_HOME

- **Tests**

To check the installation, simply run the following Scipion test:

scipion3 tests schrodingerScipion.tests.test_sitemap.TestSitemap

===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/schrodinger_devel.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/schrodinger_prod.svg
