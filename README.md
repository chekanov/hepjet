# HepJet
A library to benchmark various jet algorithms for HEP physics

All algorithms use double precision and rapidity-phi space to define distances. The directories to perform benchmarks are organized as following:

<p>
<b>fastjet</b>       - the official FastJet implementation <br>
<b>ktjet</b>         - original KtJet (C++) benchmark code. Works only in the kT mode<br>
<b>scjet_cpp</b>     - SCJet. An alternative implementation of kT-jet / anti-kT clustering in C++<br>
<b>scjet_java</b>    - Implementation of the SCJet jet algorithm in Java. <br>
<p>

FastJet and KtJet C++ codes used by this library are taken  from the orinal <a href="http://fastjet.fr/">FastJet</a> and <a href="https://ktjet.hepforge.org/">KtJet</a> web pages. 
SCjet is a light-weight implementation of the kT/anti-kT algorithms for jet validation used by the 
<a href="http://atlaswww.hep.anl.gov/hepsim/">HepSim</a> Monte Carlo database.
More details can be found in <a href="https://github.com/chekanov/hephysics">HePhysics package</a>. 



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
