import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.Collections;
import java.util.StringTokenizer;
import java.text.*;

/**
 * JetN2 is an implementation of the longitudinally-invariant kT, anti-KT and
 * Cambridge/Aachen clustering algorithms. The algorithm uses rapidity-phi for
 * the distance parameter and double values for merging. The input and output
 * for this algorithm is {@link hephysics.jet.ParticleD} class.
 * <p>
 * This class uses double values for calculations, caching and requires more
 * memory, compared to the light-weight {@link hephysics.jet.KTjet} class that
 * uses floats and pseudo-rapidity to define the distance parameter. This
 * implementation can access jet constituents.
 * <p>
 * </p>
 * This algorithm is similar to the N2 FastJet http://fastjet.fr/ implementation
 * that uses rapidity. The method uses E-scheme to combine particles (p1+p2).
 * More details is in http://arxiv.org/pdf/hep-ph/0210022v1.pdf.
 * 
 * @author Ivan Pogrebnyak <ivanp@msu.edu> and S.Chekanov (ANL)
 * 
 */
public class JetN2 {

	private int recom = 1;
	private double R;
	private final double PI2 = Math.PI * 2;
	private boolean debug = false;
	private double minpt = 0;
	private String type = "antikt";
	private ArrayList<ParticleD> jets;
	private DecimalFormat formatter = new DecimalFormat("%.12f");
	private ClusterSequence seq;

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
	 *            2: The pT-scheme (not implemented) <br>
	 *            3: The pT^2 scheme (not implemented) <br>
	 *            Currently only E-scheme is implemented.
	 * @param type
	 *            [kt,antikt,ca] <br>
	 *            clustering mode dij=min(kT_i^{2* mode},kT_j^{2* mode})). <br>
	 *            kt : means inclusive kT jet algorithm <br>
	 *            ca: means Cambridge/Aachen jet algorithm <br>
	 *            antikt: means anti-KT jet algorithm <br>
	 * @param minpt
	 *            min pT for final jets.
	 */
	public JetN2(double R, int recom, String type, double minpt) {
		this.R = R;
		this.recom = recom;
		this.debug = false;
		this.minpt = minpt;
		this.type = type.trim();

		DecimalFormat formatter1 = new DecimalFormat("#0.00");
		String rs = formatter1.format(this.R);
		System.out
				.println("JetN2: Initialization of Java jet algorithm.");
                System.out
                                .println("JetN2: Authors: S.Chekanov (ANL), I.Pogrebnyak (MSU)");

		System.out
				.println("JetN2: Inclusive mode using the E-scheme recombination and R="
						+ rs);
		if (type.equalsIgnoreCase("kt"))
			System.out.println("JetN2: Longitudinally invariant kt algorithm");
		else if (type.equalsIgnoreCase("ca"))
			System.out.println("JetN2: Cambridge/Aachen algorithm");
		else if (type.equalsIgnoreCase("antikt"))
			System.out
					.println("JetN2: Longitudinally invariant anti-kt algorithm");
		else { 
                        this.type="antikt";
			System.out
					.println("JetN2: Not correct mode:  Fallback to the inclusive kT algorithm using E-scheme and R="
							+ rs);
                     }

		if (recom != 1) {
			System.out
					.println("JetN2: Only E-scheme recombination supported! Exit.");
			System.exit(0);
		}

		seq = new ClusterSequence(this.type, R);

	}



         /**
         * Initialize calculations of the longitudinally invariant kT algorithm in
         * inclusive mode. Jet can be clustered using Cambridge/Aachen or anti-kT
         * approaches, depending on the "type" parameter. The distance parameters
         * are rapidity and phi. The E-scheme with 4-vector addition is used. 
         * 
         * @param R
         *            distance measure
         * @param type
         *            [kt,antikt,ca] <br>
         *            clustering mode dij=min(kT_i^{2* mode},kT_j^{2* mode})). <br>
         *            kt : means inclusive kT jet algorithm <br>
         *            ca: means Cambridge/Aachen jet algorithm <br>
         *            antikt: means anti-KT jet algorithm <br>
         * @param minpt
         *            min pT for final jets.
         */
        public JetN2(double R, String type, double minpt) {
             this(R, 1, type, minpt);
        }



	/**
	 * Initialize calculations of the antikT algorithm. Meaningful values are R=0.2-
	 * 1. Jets are clustered in rapidity and phi space. The E-scheme with
	 * 4-vector addition is used.
	 * 
	 * @param R
	 *            distance measure
	 * @param minpt
	 *            min pT for final jets.
	 */
	public JetN2(double R, double minpt) {
		this(R, 1, "antikt", minpt);
	}


        /**
         * Initialize calculations of the kT algorithm. Meaningful values are R=0.2-
         * 1. Jets are clustered in rapidity and phi space. The E-scheme with
         * 4-vector addition is used. Minumum pT is set to 5 GeV. 
         * 
         * @param R
         *            distance measure
         */
        public JetN2(double R) {
                this(R, 1, "kt", 5.0);
        }


	/**
	 * Run the jet algorithm using the list of particles.
	 * 
	 * @param list
	 *            list with particles
	 * @return final jets without sorting.
	 */
	public List<ParticleD> buildJets(List<ParticleD> list) {

		jets = seq.cluster(list, minpt);
                if (debug) printJets(); 

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

	

	/**
	 * Print debugging information. It shows how much time spend to make jets in
	 * ms.
	 * 
	 * @param debug
	 *            true if printing benchmark information.
	 */
	public void setDebug(boolean debug) {
		if (debug)
			System.out.println("Debug mode is ON");
		this.debug = debug;
	}


         /**
         * Return information about the jet setting. 
         * 
         * @param info 
         *            Information on the current settings. 
         */
        public String info() {


               DecimalFormat formatter1 = new DecimalFormat("#0.00");
               String rs = formatter1.format(this.R);

               String tmp="";
               tmp=tmp+"JetN2: Initialization of Java jet algorithm:\n";
               tmp=tmp+"JetN2: authors: S.Chekanov (ANL), I.Pogrebnyak (MSU)\n"; 
               tmp=tmp+"JetN2: Inclusive mode using the E-scheme recombination and R="
                                                + rs+"\n";
                if (type.equalsIgnoreCase("kt"))
                        tmp=tmp+"JetN2: Longitudinally invariant kt algorithm\n";
                else if (type.equalsIgnoreCase("ca"))
                        tmp=tmp+"JetN2: Cambridge/Aachen algorithm";
                else if (type.equalsIgnoreCase("antikt"))
                         tmp=tmp+"JetN2: Longitudinally invariant anti-kt algorithm\n";
                else {
                        this.type="antikt";
                        tmp=tmp+"JetN2: Not correct mode:  Fallback to the inclusive antikT algorithm using E-scheme and R="+ rs+"\n"; 
                     }

                if (recom != 1) {
                       
                        tmp=tmp+"JetN2: Only E-scheme recombination supported! Exit.";
                        System.exit(0);
                }

                return tmp;
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
		for (int i = 0; i < 10; ++i) {

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

			System.out.println("Number of particles="
					+ Integer.toString(list.size()));
			System.out.println("Run Nr=" + Integer.toString(i));

			long startTime = System.currentTimeMillis();
			JetN2 kt = new JetN2(0.6, "antikt", 5.0);
			kt.setDebug(false);
			kt.buildJets(list);
			kt.printJets();
			System.out.println("--->  Run time for jet creation: "
					+ Long.toString(System.currentTimeMillis() - startTime)
					+ " ms");
		}

	}

}
