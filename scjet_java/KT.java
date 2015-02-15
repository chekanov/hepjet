import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.Collections;
import java.util.StringTokenizer;
import java.text.*;

/**
 * Longitudinally-invariant kT, anti-KT and Cambridge/Aachen clustering
 * algorithms. The algorithm uses rapidity-phi for the distance parameter and
 * double values for merging. The input and output for this algorithm is
 * {@link hephysics.jet.ParticleD} class.
 * <p>
 * This class uses double values for calculations, caching and requires more
 * memory, compared to the light-weight {@link hephysics.jet.KTjet} class that
 * uses floats and pseudo-rapidity to define the distance parameter. This
 * implementation can access jet constituents.
 * <p>
 * </p>
 * This algorithm is similar to the FastJet http://fastjet.fr/ implementation
 * that uses rapidity. Use light-weight {@link hephysics.jet.KTjet} class when
 * using pseudo-rapidity and phi to define distance parameters. The method uses
 * E-scheme to combine particles (p1+p2). More details is in
 * http://arxiv.org/pdf/hep-ph/0210022v1.pdf.
 * 
 * @author S.Chekanov
 * 
 */
public class KT {

	private int recom = 1;
	private double R;
	private double R2;
	private int[] is_consider; // 0-ignore, 1-original, n>1 - merged, -1 - jet
	private double[] ktdistance1;
	private double[][] ktdistance12; // keep dij distances
	private ArrayList<ParticleD> jets;
	private final double PI2 = Math.PI * 2;
	private boolean debug = false;
	private double minpt = 0;
	private int mode = 1;
	private DecimalFormat formatter = new DecimalFormat("%.12f");

	/**
	 * Initialize calculations of the longitudinally invariant kT algorithm in
	 * inclusive mode. Jet can be clustered using Cambridge/Aachen or anti-kT
	 * approaches, depending on the "mode" parameter. The distance parameters
	 * are rapidity and phi.
	 * 
	 * @param R
	 *            distance measure
	 * @param recom
	 *            recombination scheme.<br>
	 *            1: The E-scheme Simple 4-vector addition. <br>
	 *            2: The pT-scheme. <br>
	 *            3: The pT^2 scheme. <br>
	 *            Currently only E-scheme is implemented.
	 * @param mode
	 *            clustering mode dij=min(kT_i^{2* mode},kT_j^{2* mode})). <br>
	 *            mode=1 means inclusive KT jet algorithm <br>
	 *            mode=0 means Cambridge/Aachen jet algorithm <br>
	 *            mode=-1 means anti-KT jet algorithm <br>
	 * @param minpt
	 *            min pT for final jets.
	 */
	public KT(double R, int recom, int mode, double minpt) {
		this.R = R;
		this.R2 = (R * R);
		this.recom = recom;
		this.debug = false;
		this.minpt = minpt;
		this.mode = mode;
		DecimalFormat formatter1 = new DecimalFormat("#0.00");
		String rs = formatter1.format(this.R);
		System.out
				.println("Initialization of Java jet algorithm. S.Chekanov (ANL)");
		System.out
				.println("Inclusive mode using the E-scheme recombination and R="
						+ rs);
		if (mode == 1)
			System.out.println("Longitudinally invariant kt algorithm");
		else if (mode == 0)
			System.out.println("Cambridge/Aachen algorithm");
		else if (mode == -1)
			System.out.println("Longitudinally invariant anti-kt algorithm");
		else
			System.out
					.println("Not correct mode:  Fallback to the inclusive kT algorithm using E-scheme and R="
							+ rs);

               if (recom != 1) {
                     System.out.println("Only E-E-scheme recombination supported! Exit.");
                     System.exit(0);
                }

	}

	/**
	 * Initialize calculations of the KT algorithm. Meaningful values are R=0.2-
	 * 1. Jets are clustered in rapidity and phi space.
	 * 
	 * @param R
	 *            distance measure
	 * @param recom
	 *            recombination scheme.<br>
	 *            1: The E-scheme Simple 4-vector addition. <br>
	 *            2: The pT-scheme. <br>
	 *            3: The pT^2 scheme. <br>
	 *            Currently only E-scheme is implemented.
	 * @param minpt
	 *            min pT for final jets.
	 */
	public KT(double R, int recom, double minpt) {
		this(R, recom, 1, minpt);
	}

	/**
	 * Initialize calculations of the KT algorithm. Meaningful values are R=0.2-
	 * 1. Jets are clustered in rapidity and phi space. The The E-scheme with
	 * 4-vector addition is used.
	 * 
	 * @param R
	 *            distance measure
	 * @param minpt
	 *            min pT for final jets.
	 */
	public KT(double R, double minpt) {
		this(R, 1, 1, minpt);
	}

	/**
	 * Run the jet algorithm using the list of particles.
	 * 
	 * @param list
	 *            list with particles
	 * @return final jets without sorting.
	 */
	public List<ParticleD> buildJets(List<ParticleD> list) {

		jets = new ArrayList<ParticleD>();
		int size = list.size();

		long startTime = 0;
		if (debug)
			startTime = System.currentTimeMillis();

                int i,j;

		ktdistance1 = new double[size];
		is_consider = new int[size];
		for (int m = 0; m < size; m++) {
			is_consider[m] = 1;
			ParticleD p1 = (ParticleD) list.get(m);
			ktdistance1[m] = getKtDistance1(p1);
		}

		ktdistance12 = new double[size][size];
		for (i=0; i<size - 1; i++) {
			ParticleD p1 = (ParticleD) list.get(i);
			for (j=i+1; j < size; j++) {
				ParticleD p2 = (ParticleD) list.get(j);
				ktdistance12[i][j] = getKtDistance12(p1, p2);

			}
		}

		if (debug) {
			long stopTime = System.currentTimeMillis();
			long runTime = stopTime - startTime;
			System.out
					.println("--->  Run time after making a cache of distances (ms): "
							+ runTime);
		}

		int Nstep = size;
		int iter=0;
                boolean merged=false;
                int km = -1;
                int j1 = -1;
                int j2 = -1;
                double min12 = Double.MAX_VALUE;
                double min1 = Double.MAX_VALUE;

		while (Nstep > 0) {

			min12 = Double.MAX_VALUE;
                        min1 = Double.MAX_VALUE;

                        // this is fast antiKT jet algorithm
                        // build pseudo-jet aroung particles with large pT
                        if (mode <0) { 
			// this is after reseting to a new jet
			  if (!merged) {
				for (i=0; i < size - 1; i++) {
					if (is_consider[i] <= 0)
						continue;
					for (j=i + 1; j < size; j++) {
						if (is_consider[j] <= 0)
							continue;
						if (ktdistance12[i][j] < min12) {
							min12 = ktdistance12[i][j];
							j1 = i;
							j2 = j;
						}
					}
				}
			  } else {
				// find another minimum around this jet when j1>0
				for (j=0; j < size; j++) {
					if (is_consider[j] <= 0 || j == j1)
						continue;
					if (ktdistance12[j1][j] < min12) {
						min12 = ktdistance12[j1][j];
						j1 = j1;
						j2 = j;
					}
				}

			} // end of min finding

			// find min distance to the beam
			min1 = ktdistance1[j1];
			if (ktdistance1[j2] < min1) min1 = ktdistance1[j2]; 


                        if (merged==false && Nstep==1) break;

                        }  else { 

                           // end fast antiKT
                           // start the usual kT algorithm..
                           // -----------------------------//

                     j1=0;
                     j2=0;
	             km=0;
                    // find smallest distances
                      for (i=0; i < size-1; i++) {
                        if (is_consider[i]<=0) continue;
                        for (j=i+1; j < size; j++) {
                                if (is_consider[j]<=0) continue;
                                if (ktdistance12[i][j] < min12) {
                                        min12 = ktdistance12[i][j];
                                        j1 = i;
                                        j2 = j;
                                }
                            }
                         }

                // find min distance to the beam
                     for (j = 0; j < size; j++) {
                        if (is_consider[j]<=0) continue;
                        if (ktdistance1[j] < min1) {
                                min1 = ktdistance1[j];
                                km = j;
                        }
                }


               } // end kT and CA 


			// make the decision about this particle
			merged = false;
			if (min12 < min1)  merged = true; 

			if (merged) {
				ParticleD p1 = (ParticleD) list.get(j1);
				ParticleD p2 = (ParticleD) list.get(j2);
				if (j1 != j2)
					p1.add(p2, j2); // also keeps an index
				Nstep--;
				list.set(j1, p1); // replace with p1+p2
				is_consider[j2] = 0;
				is_consider[j1] = is_consider[j1] + 1; 
				// recalculate distance for this particle
				ktdistance1[j1] = getKtDistance1(p1);
				for (i = 0; i < size; i++) {
					if (is_consider[i] <= 0 || i == j1)
						continue;
					ParticleD pp1 = (ParticleD) list.get(i);
					ktdistance12[j1][i] = getKtDistance12(p1, pp1);
				}

			}

			if (!merged) { // add this to the jet
                                if (mode>=0) j1=km; // thsi is for KT and C/A
				is_consider[j1] = -1;
				ParticleD pj = (ParticleD) list.get(j1);
				Nstep--;
				if (pj.getPt() > minpt) {
					jets.add(pj); // fill jets
				}
			}

			
			 if (debug) {
			 iter++;
                         System.out.println("## Iteration:"+Integer.toString(iter));
                         for (i=0; i< size; i++) {
                         ParticleD p1 = (ParticleD) list.get(i);
                         String mess="original";
                         if (is_consider[i]==-1) mess="!final-jet!";
                         if (is_consider[i]>1) mess="(proto-jet)";
                         if (is_consider[i]==0) mess="(removed)";
                        System.out.println( Integer.toString(i)+"  E="+Double.toString(p1.e())+" "+mess);
                 }
         }

			
			// end loop
		}

		
	if (debug) {

        System.out.println("Final Nr of iterations="+Integer.toString(iter));

		
        // attempt to deal with unmeargable particle
        int ins=-1;
        for (i=0; i < size; i++)
        if (is_consider[i]==1) {ins=i;};

        if (ins>-1) {
                ParticleD p2 = list.get(ins);
                if (debug) System.out.println("Unmerged particle id="+Integer.toString(ins));
                min12 = Double.MAX_VALUE;
                for (j=0; j < jets.size(); j++) {
                        ParticleD lp = jets.get(j);
                        double d=getDistance(p2, lp);
                        if (d<min12) { j1=j; min12 =d; };
                }
                if (debug)  System.out.println("Distance R to closest jet="+Double.toString(min12));
                if (min12<R) {
                        if (debug)  System.out.println(" --> Particle merged");
                        ParticleD lp = jets.get(j1);
                        lp.add(p2,j1);
                        is_consider[ins] = 0;
                }
        }

                // sanity test. All particles were merged?
                int nn=0; ins=-1;
                for (i=0; i<size; i++)
                if (is_consider[i]==1) {nn++; ins=i;};
                if (nn != 0)   System.out.println( "--> WARNING: particle with ID="+ Integer.toString(ins)+" unmerged");


			long stopTime2 = System.currentTimeMillis();
			long runTime = stopTime2 - startTime;
			System.out.println("  --> Final time for calculation (ms): "
					+ runTime);
			System.out.println("  --> Nr of jets : " + jets.size());

	} // end debug mode

		is_consider = null;
		ktdistance12 = null;
		ktdistance1 = null;
		return jets;

	}

	/**
	 * Get jets after sorting in jet pT. Run buildJets before calling this
	 * method.
	 * 
	 * @return list with sorted jets
	 */

	public ArrayList<ParticleD> getJetsSorted() {

		Collections.sort(jets);

		return jets;
	}

	/**
	 * Print the kT jets for debugging.
	 */
	public void printJets() {

		ArrayList<ParticleD> sjets = getJetsSorted();

		System.out.println("# Nr of jets=" + Integer.toString(sjets.size()));
		System.out.format("%5s %14s %14s %14s %7s\n", "jet #", "rapidity",
				"phi", "pt", " const");
		for (int i = 0; i < sjets.size(); i++) {
			ParticleD lp = sjets.get(i);
			double phi = lp.phi();
			List con = lp.getConstituentsList();
			if (phi < 0)
				phi = PI2 + phi;
			String s1 = String.format("%15.8f", lp.getRapidity());
			String s2 = String.format("%15.8f", phi);
			String s3 = String.format("%15.8f", lp.getPt());
			String nc = Integer.toString(con.size());
			System.out.format("%5s%15s%15s%15s%7s\n", Integer.toString(i), s1,
					s2, s3, nc);
		}

	}

	/**
	 * Print the kT jets for debugging to a string.
	 * 
	 * @return String representing a jet
	 */
	public String toString() {
		ArrayList<ParticleD> sjets = getJetsSorted();
		String tmp = "# Nr of jets=" + Integer.toString(sjets.size()) + "\n";
		for (int i = 0; i < sjets.size(); i++) {
			ParticleD lp = (ParticleD) sjets.get(i);
			List con = lp.getConstituentsList();
			String spx = formatter.format(lp.getRapidity());
			String spy = formatter.format(lp.getPhi());
			String spz = formatter.format(lp.getPt());
			tmp = tmp + "n=" + Integer.toString(i) + " y=" + spx + " phi="
					+ spy + " pt=" + spz + " const="
					+ Integer.toString(con.size()) + "\n";
		}

		return tmp;
	}

	private double phiAngle(double phi) {
		if (phi > PI2)
			phi -= (PI2);
		if (phi < -PI2)
			phi += (PI2);
		return phi;
	}

	/**
	 * Calculate delta R distance.
	 * 
	 * @param a
	 *            input particle
	 * @param b
	 *            input particle
	 * @param p
	 *            power parameter
	 * @return Kt distance
	 */
	public double getKtDistance12(ParticleD a, ParticleD b) {
		double rsq, esq, deltaEta, deltaPhi;
		deltaEta = a.getRapidity() - b.getRapidity();
		double phi1 = a.getPhi();
		double phi2 = b.getPhi();
		deltaPhi = phi2 - phi1;
		if (deltaPhi > Math.PI)
			deltaPhi = PI2 - deltaPhi;
		if (deltaPhi < -Math.PI)
			deltaPhi = PI2 + deltaPhi;
		rsq = (deltaEta * deltaEta + deltaPhi * deltaPhi);
		esq = 0;
		if (mode == 1)
			esq = Math.min(a.getPt2(), b.getPt2()); // kT
		else if (mode == 0)
			esq = 1.0;  // C-A
		else if (mode == -1)
			esq = Math.min(1.0 / a.getPt2(), 1.0 / b.getPt2()); // antiKT
		else
			esq = Math.min(a.getPt2(), b.getPt2()); // kT

		return (esq * rsq / R2);
	}


         /**
         * Calculate  R distance in y-phi.
         * 
         * @param a
         *            input particle
         * @param b
         *            input particle
         * @return y-phi distance
         */
        public double getDistance(ParticleD a, ParticleD b) {
                double rsq, deltaEta, deltaPhi;
                deltaEta = a.getRapidity() - b.getRapidity();
                double phi1 = a.getPhi();
                double phi2 = b.getPhi();
                deltaPhi = phi2 - phi1;
                if (deltaPhi > Math.PI)
                        deltaPhi = PI2 - deltaPhi;
                if (deltaPhi < -Math.PI)
                        deltaPhi = PI2 + deltaPhi;
                rsq = (deltaEta * deltaEta + deltaPhi * deltaPhi);
                return Math.sqrt(rsq);
        }

	/**
	 * This is the KT distance to the beam (assuming Z=Y=0). The distance
	 * measure depends on the mode parameter.
	 * 
	 * @param a
	 *            particle
	 * @return kT distance
	 */
	public double getKtDistance1(ParticleD a) {
		if (mode == 1)
			return a.getPt2();
		else if (mode == 0)
			return 1.0;
		else if (mode == -1)
			return (1.0 / a.getPt2());
		return a.getPt2();
	}

	/**
	 * Print debugging information. It shows how much time spend to make jets in
	 * ms.
	 * 
	 * @param debug
	 *            true if printing benchmark information.
	 */
	public void setDebug(boolean debug) {
                if (debug) System.out.println("Debug mode is ON");
		this.debug = debug;
	}

	/**
	 * Main class for testing.
	 * 
	 * @param args
	 */

	public static void main(String[] args) {

		String data = "";
		if (args.length > 0) {
			data = args[0];
		} else {
			System.out.println("No input file with particles! Exit!");
			System.exit(1);

		}

		// for correct benchmark with C++ (after just-in-time compiler)
		for (int i = 0; i < 3; i++) {

			List<ParticleD> list = new ArrayList<ParticleD>();
			try {
				File file = new File(data);
				FileReader fileReader = new FileReader(file);
				BufferedReader bufferedReader = new BufferedReader(fileReader);
				String line;
				while ((line = bufferedReader.readLine()) != null) {

					StringTokenizer st = new StringTokenizer(line);
					int j = 0;
					double[] mom = new double[4];
					while (st.hasMoreElements()) {
						Double d = Double.parseDouble(st.nextElement()
								.toString());
						mom[j] = d;
						j++;
					}

					// px,py,pz,e
					ParticleD pp = new ParticleD(mom[0], mom[1], mom[2], mom[3]);
					list.add(pp);

				}
				fileReader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

			System.out.println("Number of particles="+Integer.toString(list.size())); 
			System.out.println("Run Nr=" + Integer.toString(i));

                        long startTime = System.currentTimeMillis();
			KT kt = new KT(0.6, 1, -1, 5.0);
			kt.setDebug(false);
			kt.buildJets(list);
			kt.printJets();
                        System.out.println("--->  Run time for jet creation: "+Long.toString(System.currentTimeMillis()-startTime)+ " ms"); 
		}

	}

}
