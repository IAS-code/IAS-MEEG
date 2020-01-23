.. IAS documentation master file, created by
   sphinx-quickstart on Wed Jan 30 17:57:08 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

IAS algorithm
=============

The Iterative Alternating Sequential (IAS) [Ref1]_ [Ref2]_ algorithm is based on an iterative scheme that alternatively updates the dipole
moments Q by solving a linear least squares problem using a priorconditioned CGLS algorithm with sutable 
stopping condition and updating the hyperparameter theta by an explicit formula.
   
An example of application of IAS to a real dataset can be found in :ref:`example`.   


Download
========


.. code-block:: bash

    git clone https://github.com/IAS-code/IAS-MEG.git 
    

.. [Ref1] D. Calvetti, A. Pascarella, F. Pitolli, E. Somersalo, B. Vantaggi
          *A hierarchical Krylov-Bayes iterative inverse solver for MEG with physiological preconditioning*, Inverse Problems, 31 (12), 125005, (2015)
 
.. [Ref2] D. Calvetti, A. Pascarella, F. Pitolli, E. Somersalo, B. Vantaggi,
          *Brain Activity Mapping from MEG Data via a Hierarchical Bayesian Algorithm with Automatic Depth Weighting*, Brain topography, 1-31, (2018) 

   
.. toctree::
   :maxdepth: 1

   example
   api
   



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
