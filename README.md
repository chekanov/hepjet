# HepJet
A library to benchmark various jet algorithms for HEP physics

All algorithms use double precision and rapidity-phi space to define distances. The directories to perform benchmarks are organized as following:

<p>
<b>fastjet</b>       - the official FastJet implementation <br>
<b>ktjet</b>         - a modified KtJet (C++) benchmark code. Works in the anti-KT  mode<br>
<b>nlojet</b>        - implementation of kT-jets from NLOjet++ <br>
<b>scjet_cpp</b>     - SCJet. An alternative implementation of kT-jet / anti-kT clustering in C++<br>
<b>scjet_java</b>    - implementation of the SCJet jet algorithm in Java. <br>
<b>benchmark</b>    -  compare FastJet and SCJet implementations using the same input.<br>
<p>
All algorithms run the anti-kT jet algorithm. One can also run the standard kT and Cambridge/Aachen jet algorithms. 

<p>

FastJet and KtJet C++ codes used by this library are taken  from the orinal <a href="http://fastjet.fr/">FastJet</a> and <a href="https://ktjet.hepforge.org/">KtJet</a> web pages. 
SCjet is a light-weight implementation of the anti-kT algorithms for jet validation used by the 
<a href="http://atlaswww.hep.anl.gov/hepsim/">HepSim</a> Monte Carlo database.
More details can be found in <a href="https://github.com/chekanov/hephysics">HePhysics package</a>. 
The Java implementation of SCJet is available from the <a href="http://jwork.org/scavis/">SCaVis data-analysis</a> program. 

<h2>Benchmark results</h2>

The benchmarks have been done on Xeon(R) CPU E5520 @ 2.27GHz
using the input file with 1000 particles (data/single-event.dat) for a 100 pp collision event
(ttbar, pT>8 TeV). The processing time:

 <ul>
  <li>fastjet    - 3 msec </li>
  <li>scjet_cpp -  105 msec (fast N^2 mode) or 990 msec (the standard N^3 mode) </li>   
  <li>scjet_java - 103 msec (after multiple runs, fast mode) </li>
  <li>nlojet    - ? (failed). Used by NLOjet++</li>
  <li>ktjet     - 370 msec</li>
  <li>pgjet_cpp - 4756 ms msec (no cache, using std::set)</li>
</ul> 

In summary: <b>fastjet</b> is  about a factor 30 faster than <b>scjet_cpp</b> (seeded).
The Java implementation (<b>scjet_java</b>) is as fast as the C++ version when using more than one run           
over events (first run is a factor 4 slower than for the C++ version due to the JIT compilation).
Other similar algorithms are slower. When using the kT mode, the fast mode is disabled
and <b>scjet_java</b> has the same speed as <b>ktjet</b>. 

<p>
</p>

There are some differences in the output jets between different implementations. 
There is no difference between <b>scjet_cpp</b> and <b>fastjet</b> implementations.
For the seeded option,  
differences can be at the level of 2% for transverse momentum of very soft jets.
The difference between <b>nlojet</b> and  <b>fastjet</b> is also at the level of a few percents. 
No difference is found between <b>fastjet</b> and  <b>ktjet</b> implementations. 

<p>


Read "README" files in each directory to see how to run each benchmark. 
Makefile files are included inside each directory.

Authors: 
--------
S.Chekanov (ANL)
I.Pogrebnyak (MSU/ANL)

Feb. 2015

