Frequently-Asked Questions
==========================

wien2k: FERMI ERROR when running `x lapw2 -almd -band`
------------------------------------------------------
In some versions of Wien2k, there is a problem in running `x lapw2 -almd -band`.

A hack solution is as follows:
1) `x lapw1 -band`
2) edit in2 file: replace 'TOT' with 'QTL', 'TETRA' with 'ROOT'
3) `x lapw2 -almd -band`
4) `dmftproj -band` (add the fermi energy to file, it can be found by running `grep :FER *.scf`)

How do I plot the output of `spaghettis`?
-----------------------------------------

In python, you can do the following for example. You should pass the name of
the file written out by the spaghettis function.  Of course, you should change
the parameters as desired.

.. literalinclude:: plotting_spaghettis.py

Why is my calculation not working?
----------------------------------

Are you running in the right shell?
