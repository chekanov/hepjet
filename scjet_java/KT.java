import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.Collections;
import java.util.StringTokenizer;
import java.text.*;

/**
 * Longitudinally-invariant kT, anti-KT and  Cambridge/Aachen clustering algorithms. 
 * The algorithm uses rapidity-phi for the distance parameter and double values for merging. The input and output for this algorithm  is  {@link hephysics.jet.ParticleD} class. <p> 
 * This class uses double values for calculations, caching and requires more memory, compared to the light-weight {@link hephysics.jet.KTjet} class that uses floats and pseudo-rapidity to define the distance parameter. This implementation can access jet constituents.  <p></p> 
 * This algorithm is similar to the FastJet http://fastjet.fr/ implementation that uses rapidity.  Use light-weight {@link hephysics.jet.KTjet} class when using pseudo-rapidity and phi to define distance parameters. The method uses E-scheme to combine particles (p1+p2). More details is in http://arxiv.org/pdf/hep-ph/0210022v1.pdf. 
 * 
 * @author S.Chekanov
 * 
 */
public class KT {

	private int recom = 1;
	private static double R;
	private static double R2;
	private boolean[] is_consider;
	private double[] ktdistance1;
	private double[][] ktdistance12;
	private ArrayList<ParticleD> jets;
	static private final double PI2 = Math.PI * 2;
	private boolean debug=false;
	private double minpt=0;
	private static int mode=1;
	private DecimalFormat formatter = new DecimalFormat("#.#####");


	/**
	 * Initialize calculations of the longitudinally invariant kT algorithm in inclusive mode. 
	        * Jet can be clustered using Cambridge/Aachen or anti-kT approaches, depending on the "mode" parameter. 
	        * The distance parameters are rapidity and phi. 
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
	        *          clustering mode dij=min(kT_i^{2* mode},kT_j^{2* mode})). <br>
	        *          mode=1 means inclusive KT jet algorithm <br> 
	        *          mode=0 means Cambridge/Aachen jet algorithm <br> 
	        *          mode=-1 means anti-KT jet algorithm <br> 
	        * @param minpt
	        *            min pT for final jets.
	 */
	public KT(double R, int recom, int mode, double minpt) {
		this.R = R;
		this.R2 = (R * R);
		this.recom = recom;
		this.debug=false;
		this.minpt=minpt;
		this.mode=mode;
		DecimalFormat formatter = new DecimalFormat("#0.00");
		String rs=formatter.format(this.R);
		System.out.println("Initialization of Java jet algorithm. S.Chekanov (ANL)");
		System.out.println("Inclusive mode using the E-scheme recombination and R="+rs);
		if (mode==1) System.out.println("Longitudinally invariant kt algorithm");
		else if (mode==0) System.out.println("Cambridge/Aachen algorithm" );
		else if (mode==-1) System.out.println("Longitudinally invariant anti-kt algorithm");
		else  System.out.println("Not correct mode:  Fallback to the inclusive kT algorithm using E-scheme and R="+rs);
	}


	/** Initialize calculations of the KT algorithm. Meaningful values are R=0.2- 1.
	* Jets are clustered in rapidity and phi space. 
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
		this(R,recom,1, minpt);
	}


	/** Initialize calculations of the KT algorithm. Meaningful values are R=0.2- 1.
	* Jets are clustered in rapidity and phi space. The The E-scheme with 4-vector addition is used. 
	* 
	* @param R
	*            distance measure
	* @param minpt
	*            min pT for final jets.
	*/
	public KT(double R, double minpt) {
		this(R,1, 1, minpt);
	}




	/**
	  * Run the jet algorithm using the list of particles. 
	  * 
	  * @param list
	  *            list with particles
	  * @return final jets without sorting.
	  */
	public  List<ParticleD>  buildJets(List<ParticleD> list) {

		jets = new ArrayList<ParticleD>();
		int size = list.size();

		long startTime = 0;
		if (debug) startTime=System.currentTimeMillis();


		ktdistance1 = new double[size];
		is_consider = new boolean[size];
		for (int m = 0; m < size; m++) {
			is_consider[m] = true;
			ParticleD p1 = (ParticleD) list.get(m);
			ktdistance1[m] = getKtDistance1(p1);
		}


		ktdistance12 = new double[size][size];
		for (int i = 0; i < size - 1; i++) {
			ParticleD p1 = (ParticleD)list.get(i);
			for (int j = i + 1; j < size; j++) {
				ParticleD p2 = (ParticleD) list.get(j);
				ktdistance12[i][j] = getKtDistance12(p1, p2);
				//System.out.println(ktdistance12[i][j]);
			}
		}

		if (debug) {
			long stopTime = System.currentTimeMillis();
			long runTime = stopTime - startTime;
			System.out.println("--->  Run time after making a cache of distances (ms): " + runTime);
		}



		int Nstep = size;
		// start loop over all objects
		int ntot=size;
		while (Nstep > 0) {

			int j1 = 0;
			int j2 = 0;
			double min12 = Double.MAX_VALUE;
			//
			// find smallest d12.
			//
			for (int i = 0; i < size-1; i++) {
				if (is_consider[i] == false)
					continue;
				for (int j = i+1; j < size; j++) {
					if (is_consider[j] == false)
						continue;
					if (ktdistance12[i][j] < min12) {
						min12 = ktdistance12[i][j];
						j1 = i;
						j2 = j;
					}
				}
			}

			// find min distance to the beam
			int km = 0;
			double min1 = Double.MAX_VALUE;
			for (int i = 0; i < size; i++) {
				if (is_consider[i] == false) continue;
				if (ktdistance1[i] < min1) {
					min1 = ktdistance1[i];
					km = i;
				}
			}

			// remove from consideration pair
			if (min12>min1) {
				is_consider[km] = false;
				Nstep--;
				ParticleD pj = (ParticleD) list.get(km);
				if (pj.getEt() > minpt) {
					jets.add(pj); // fill jets
					// System.out.println(pj.getPt() + " " + minpt);
				}
			} else { // when min1min12< combine

				ParticleD p1 = (ParticleD) list.get(j1);
				ParticleD p2 = (ParticleD) list.get(j2);
				if (j1 != j2) p1.add(p2,j2); // also keeps an index
				Nstep--;
				list.set(j1, p1); // replace with p1+p2
				// list.remove(j2); // remove softest
				is_consider[j2] = false; // remove softest
				// recalculate distance for this particle
				ktdistance1[j1]  = getKtDistance1(p1);
				for (int i = 0; i < size; i++) {
					if (is_consider[i] == false)  continue;
					ParticleD pp1 = (ParticleD) list.get(i);
					ktdistance12[j1][i] = getKtDistance12(p1, pp1);
				}

			}



			// end loop
		}


		if (debug) {
			long stopTime2 = System.currentTimeMillis();
			long runTime = stopTime2 - startTime;
			System.out.println("  --> Final time for calculation (ms): " + runTime);
			System.out.println("  --> Nr of jets : " + jets.size());

		}

		is_consider=null;
		ktdistance12=null;
		ktdistance1=null;
		return jets;

	}




	/**
	* Get jets after  sorting in jet pT. Run  buildJets before calling this method. 
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

		ArrayList<ParticleD> sjets=getJetsSorted();

		System.out.println("# Nr of jets=" + Integer.toString(sjets.size()));
		for (int i = 0; i < sjets.size(); i++) {
			ParticleD lp = sjets.get(i);
			double phi = lp.phi();
			List con=lp.getConstituentsList();
			if (phi<0) phi=PI2+phi;
			String s1=formatter.format(lp.getRapidity());
			String s2=formatter.format(phi);
			String s3=formatter.format(lp.getEt());
			System.out.println("n=" + Integer.toString(i) + " y="
			                   + s1 + " phi="
			                   + s2 + " pt="
			                   + s3 + " const=" + Integer.toString(con.size()));

		}

	}


	/**
	* Print the kT jets for debugging to a string.
	* @return String representing a jet 
	*/
	public String toString() {
		ArrayList<ParticleD> sjets=getJetsSorted();
		String tmp="# Nr of jets=" + Integer.toString(sjets.size())+"\n";
		for (int i = 0; i < sjets.size(); i++) {
			ParticleD lp = (ParticleD)sjets.get(i);
			List con=lp.getConstituentsList();
			String spx=formatter.format(lp.getRapidity());
			String spy=formatter.format(lp.getPhi());
			String spz=formatter.format(lp.getEt());
			tmp=tmp+"n="+Integer.toString(i)+" y="+spx+" phi="+spy+" pt="+spz+" const=" + Integer.toString(con.size())+"\n";
		}

		return tmp;
	}




	private static double phiAngle(double phi) {
		if (phi>PI2) phi -= (PI2);
		if (phi<-PI2) phi += (PI2);
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
	public static double getKtDistance12(ParticleD a, ParticleD b) {
		double rsq, esq, deltaEta, deltaPhi;
		deltaEta = a.getRapidity() - b.getRapidity();
		double phi1 = a.getPhi();
		double phi2 = b.getPhi();
		deltaPhi =  phi1-phi2;
		rsq = (deltaEta*deltaEta + deltaPhi*deltaPhi);
		esq = 0;
		if      (mode==1) esq=Math.min(a.getEt2(), b.getEt2());         // kT
		else if (mode==0) esq=Math.min(a.getEt(), b.getEt());           // C-A
		else if (mode==-1) esq=Math.min(1.0/a.getEt2(), 1.0/b.getEt2()); // antiKT
		else esq=Math.min(a.getEt2(), b.getEt2());         // kT

		return (esq * rsq/R2);
	}

	/**
	* This is the KT distance to the beam (assuming Z=Y=0). 
	* The distance measure depends on the mode parameter. 
	* @param a
	*            particle
	* @return kT distance
	*/
	public static  double getKtDistance1(ParticleD a) {
		if (mode==1) return a.getEt2();
		else if (mode==0) return a.getEt();
		else if (mode==-1) return (1.0/a.getEt2());
		return a.getEt2();
	}



	/**
	 * Print debugging information. It shows how much time spend to make jets in ms.
	 * @param debug true if printing benchmark information. 
	 */
	public void setDebug(boolean debug){

		this.debug=debug;
	}

	/**
	 * Main class for testing.
	 * 
	 * @param args
	 */

	public static void main(String[] args) {


                String data="";
                if (args.length > 0) {
                   data=args[0];
                } else {
                   System.out.println("No input file with particles! Exit!");
                   System.exit(1); 
 
                   }

		List<ParticleD> list = new ArrayList<ParticleD>();
		try {
			File file = new File(data);
			System.out
			.println("Reading test file with jets. Number of particles:");
			FileReader fileReader = new FileReader(file);
			BufferedReader bufferedReader = new BufferedReader(fileReader);
			String line;
			while ((line = bufferedReader.readLine()) != null) {

				StringTokenizer st = new StringTokenizer(line);
				int j = 0;
				double[] mom = new double[4];
				while (st.hasMoreElements()) {
					Double d = Double.parseDouble(st.nextElement().toString());
					mom[j] = d;
					j++;
				}

				// px,py,pz,e
				ParticleD  pp = new  ParticleD(mom[0], mom[1],mom[2],mom[3]);
				//System.out.println(pp.perp());
				//System.out.println(pp.getEt());
				list.add(pp);

			}
			fileReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// System.out.println(list.size());
		KT KT = new KT(0.6, 1, -1, 5.0);
		KT.setDebug(true);
		KT.buildJets(list);
		KT.printJets();


		//		       Ran Longitudinally invariant kt algorithm with R = 0.6 and E scheme recombination
		//             fastjet implementation
		//				jet #        rapidity             phi              pt
		//				    0     -0.86706550      2.90526958    983.06499827
		//				    1      0.22086793      6.02950080    898.49716062
		//				    2     -1.17103708      6.07159548     69.71050689
		//				    3      0.36257181      0.54763129     16.65531147
		//				    4     -2.46859216      1.03398974      7.98755167
		//				    5     -1.63742088      4.01894081      7.60023410


	}

}
