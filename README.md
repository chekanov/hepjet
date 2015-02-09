# HepJet
A library to benchmak verious jet algorithms for HEP physics

All algorithms use double precision and rapidity-phi space to define distances. The directories to performe benchmarks are organized as following:

<p>
<b>fastjet</b>       - the official FastJet implementation <br>
<b>ktjet</b>         - original KtJet++ benchmark code. Works only in the kT mode<br>
<b>scjet_cpp</b>     - an alternative implementation of kT-jet / anti-kT clustering <br>
<b>scjet_java</b>    - Java implementation of the previous c++ algorithm. <br>
<p>


FastJet is typically a factor 40 faster than scjet_cpp and factor 160 faster than the Java implementation of the simple algorithm. There are however some differences in the output jets between different implementations.  No difference is detected between the FastJet and KtJet implementations using the kT mode (R=0.6). 

<p>


Read "README" file to see how to run each case. Makefiles are included inside each directory.

S.Chekanov (ANL)
