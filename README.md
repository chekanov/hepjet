# HepJet
A small library to benchmak jet algorithms for HEP physics

All algorithms use double precision and rapidity-phi space to define distances. The directories are organized as following:

<p>
<b>fastjet</b>       - the official FastJet implementation <br>
<b>scjet_cpp</b>     - a simple (slow) implementation of KT-jet algorithms <br>
<b>scjet_java</b>   - Java implementation of the previous simple algorithm. <br>
<b>ktjet</b>         - KtJet++ benchmark code. Works only for the kT mode<br>
<p>


FastJet is typically a factor 40 faster than scjet_cpp and factor 160 faster than the Java implementation of the simple algorithm. There are however some differences in the output jets between different implementations.  No difference is detected between fastjet and KtJet when using the kT mode. 

<p>


Read "README" file to see how to run each case. Makefiles are included inside each directory.

S.Chekanov (ANL)
