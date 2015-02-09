# HepJet
A library to benchmark various jet algorithms for HEP physics

All algorithms use double precision and rapidity-phi space to define distances. The directories to perform benchmarks are organized as following:

<p>
<b>fastjet</b>       - the official FastJet implementation <br>
<b>ktjet</b>         - original KtJet (C++) benchmark code. Works only in the kT mode<br>
<b>scjet_cpp</b>     - an alternative implementation of kT-jet / anti-kT clustering in C++<br>
<b>scjet_java</b>    - Java implementation of the previous c++ algorithm. <br>
<p>


<h2>Benchmark results</h2>
FastJet is  about a factor 40 faster than SCJet (<b>scjet_cpp</b>) and a factor 160 faster than the Java implementation (<b>scjet_java</b>) of the SCJet.  KtJet and SCJet in C++ have similar runtime performance.  

<p>
</p>

However, there are some differences in the output jets between different implementations. 
The difference between SCJet and FastJet implementations 
is at the level of 2-10% for subleading jets. 
No difference was found between the FastJet and KtJet C++ implementations using the kT mode (R=0.6). 

<p>


Read "README" files in each directory to see how to run each benchmark. 
Makefile files are included inside each directory.

S.Chekanov (ANL)
