# hepjet
A small library to benchmak jet algorithms for HEP physics

All algorithms use double precision and eta-phi space to define distances.
<<p>

fastjet       - the official FastJet implementation <br>
scjet_cpp     - a simple (slow) implementation of KT-jet algorithms <br>
scijet_java   - Java implementation of the previous simple algorithm. <br>

<p>

FastJet is typically a factor 40 faster than scjet_cpp,
and factor 160 faster than the Java implementation of the simple algorithm. There are however some
differences in the output jets.

<p>


Read "README" file to see how to run each case.

S.Chekanov (ANL)

