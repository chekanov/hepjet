# hepjet
A small library to benchmak jet algorithms for HEP physics

All algorithms use double precision and eta-phi space to define distances.

fastjet       - the official FastJet implementation
scjet_cpp     - a simple (slow) implementation of KT-jet algorithms
scijet_java   - Java implementation of the previous simple algorithm.

FastJet is typically a factor 40 faster than scjet_cpp,
and factor 160 faster than the Java implementation of the simple algorithm. There are however some
differences in the output jets.

Read "README" file to see how to run each case.

S.Chekanov

