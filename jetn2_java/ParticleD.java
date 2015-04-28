import java.io.Serializable;
import java.text.*;
import java.util.*;

/**
 * A class representing a jet or particle with pre-computed px,py,pz,e (double
 * precision). It uses double types to keep information on 4-meomenta. The
 * merging is done in the 4-momentum space. The class has a minimum dynamic
 * computation to minimize CPU usage.
 * 
 * @author S.Chekanov
 * 
 */
public class ParticleD implements Comparable<ParticleD>, Serializable {

	private static final long serialVersionUID = 1L;

	private double px, py, pz;
	private double energy;
	private final double PI2 = Math.PI * 2;
	private double phi, pt2, rapidity;
	private DecimalFormat formatter = new DecimalFormat("0.###E0");
	private ArrayList<Integer> consts;

	/**
	 * Initialize pseudoparticle.
	 * 
	 */
	public ParticleD() {
		consts = new ArrayList<Integer>();
	}

	/**
	 * Initialize from 4-momenta.
	 * 
	 * @param px
	 * @param py
	 * @param pz
	 * @param energy
	 */
	public ParticleD(double px, double py, double pz, double energy) {
		this.px = px;
		this.py = py;
		this.pz = pz;
		this.energy = energy;
		cachePhiRapidity();
		consts = new ArrayList<Integer>();
		consts.add(new Integer(-1));
	}

	/**
	 * Initialize from 4-momenta.
	 * 
	 * @param px
	 * @param py
	 * @param pz
	 * @param energy
	 */
	public void setPxPyPzE(double px, double py, double pz, double energy) {
		this.px = px;
		this.py = py;
		this.pz = pz;
		this.energy = energy;
		cachePhiRapidity();
	}

	/**
	 * Compute pseudorapidity.
	 * 
	 * @return  pseudorapidity. 
	 */
	public double eta() {
		if (px == 0.0 && py == 0.0)
			return -999;
		if (pz == 0.0)
			return -999;
		double pt2 = px * px + py * py;
		double theta = Math.atan2(Math.sqrt(pt2), pz);
		if (theta < 0)
			theta += Math.PI;
		return -1 * Math.log(Math.tan(theta / 2));
	}

	/**
	 * Compute transverse energy squared.
         * @return transverse energy squared.
	 */
	public double et2() {
		double pt2x = perp2();
		double et2 = pt2x == 0 ? 0 : e() * e() * pt2 / (pt2x + pz() * pz());
		return et2;
	}

	/**
	 * Compute transverse energy.
	 * 
	 * @return transverse energy. 
	 */
	public double et() {
		double etet = et2();
		return e() < 0.0 ? -Math.sqrt(etet) : Math.sqrt(etet);
	}

        /**
         * Compute mass.
         * 
         * @return mass 
         */
        public double mass() {
                double m=energy*energy-px*px-py*py-pz*pz;
                if (m>=0) return Math.sqrt(m);
                return -1; 
        }

        /**
         * Compute mass.
         * 
         * @return mass 
         */
        public double m() {
                return mass();
        }

	/**
	 * Compute rapidity. 0.5*log( (m+z)/(m-z) );
	 * 
	 * @return rapidity 
	 */
	public double rapidity() {
		rapidity = -10e10;
		if (energy > pz())
			rapidity = 0.5 * Math.log((energy + pz) / (energy - pz));
		return rapidity;
	}

	/**
	 * Compute magnitude sqrt(px**2+py**2+pz**2)
	 * 
	 * @return mag 
	 */
	public double mag() {
		return Math.sqrt(px * px + py * py + pz * pz);
	}

	/**
	 * Compute pT**2. 
	 * 
	 * @return pt**2 
	 */
	public double perp2() {
		return (px * px + py * py);
	}

	/**
	 * Compute transverse momentum (pT). 
	 * 
	 * @return Transverse momentum (pt) 
	 */
	public double perp() {
		return Math.sqrt(perp2());
	}


         /**
         * Set energy.
         * @return Transverse momentum (pt) 
         */
	public void setEnergy(double energy) {
		this.energy = energy;

	}

	/**
	 * Comparator. Use pT2  for comparison (in increasing order)
	 * 
	 * @param o particle 
	 * @return
	 */
	public int compareTo(ParticleD o) {
		if (perp2() < o.perp2())
			return 1;
		if (perp2() > o.perp2())
			return -1;
		return 0;
	}

	/**
	 * Copy a particle.
	 * 
	 * @param o
	 * @return a copy of particle
	 */
	public ParticleD copy(ParticleD o) {
		ParticleD tmp = new ParticleD(px, py, pz, energy);
		tmp.setRapidity(rapidity);
		tmp.setPhi(phi);
		tmp.setPt2(pt2);
		return tmp;
	}

	/**
	 * Convert a particle to a string.
	 * 
	 * @return a string with the particle
	 */
	public String toString() {
		String spx = formatter.format(px);
		String spy = formatter.format(py);
		String spz = formatter.format(pz);
		String se = formatter.format(energy);
		String srap = formatter.format(rapidity);
		String sphi = formatter.format(phi);
		String set = formatter.format(Math.sqrt(pt2));
		return "px=" + spx + " py=" + spy + " pz=" + spz + " e=" + se + " y="
				+ srap + " phi=" + sphi + " pt=" + set;
	}

	/**
	 * Get px.
	 * 
	 * @return px component
	 */
	public double px() {
		return px;
	}

	/**
	 * Get py.
	 * 
	 * @return py component
	 */
	public double py() {
		return py;
	}

	/**
	 * Get pz.
	 * 
	 * @return pz component
	 */
	public double pz() {
		return pz;
	}

	/**
	 * Get cached rapidity
	 * 
	 * @return rapidity
	 */
	public double getRapidity() {
		return rapidity;
	}

	/**
	 * Get cached perp**2.
	 * 
	 * @return perp2.
	 */
	public double getPt2() {
		return pt2;
	}

	/**
	 * Get cached Pt.
	 * 
	 * @return et transverse energy
	 */
	public double getPt() {
		return Math.sqrt(pt2);
	}

	/**
	 * Get cached phi
	 * 
	 * @return cached phi
	 */
	public double getPhi() {
		return phi;
	}

	/**
	 * Get indexes of constituents, iff filled with the "add" method
	 * 
	 * @return list of constituents.
	 */
	public List<Integer> getConstituentsList() {
		return consts;
	}

	public void setConstituents(ArrayList<Integer> consts) {
		this.consts = consts;
	}

	public void addConstituent(int i) {
		consts.add(i);
	}

	/**
	 * Get indexes of constituents, iff filled with the "add" method
	 * 
	 * @return array of constituents.
	 */
	public int[] getConstituents() {
		int s = consts.size();
		int[] intArray = new int[s];
		for (int i = 0; i < s; i++) {
			intArray[i] = consts.get(i).intValue();
		}

		return intArray;
	}

	/**
	 * Compute Phi
	 * 
	 * @return
	 */
	public double phi() {
		if (px == 0)
			return 0.0;
		phi = Math.atan2(py, px);
		// if (phi<0) phi = PI2+phi;
		return phi;
	}

	/**
	 * Get energy.
	 * 
	 * @return energy component
	 */
	public double e() {
		return energy;
	}

	public void setPx(double px) {
		this.px = px;
	}

	public void setPy(double py) {
		this.py = py;
	}

	public void setPz(double pz) {
		this.pz = pz;
	}

	public void setRapidity(double rapidity) {
		this.rapidity = rapidity;
	}

	public void setPhi(double phi) {
		this.phi = phi;
	}

	public void setPt2(double pt2) {
		this.pt2 = pt2;
	}

	/**
	 * The method precomputers Phi, Rapidity and Pt2 and store them. Such
	 * caching makes faster computations. Use getRapidity(), getPhi(), getPt2()
	 * methods to return such values.
	 * 
	 */
	public void cachePhiRapidity() {
		rapidity = rapidity();
		phi = phi();
		pt2 = perp2();
	}

	/**
	 * Add to this particle another particle. The method also precomputes
	 * rapidity, phi and pt2 for fast retrival using getPhi, getRapidity
	 * methods.
	 * 
	 * @param a
	 */
	public void add(ParticleD a) {
		px = a.px() + px;
		py = a.py() + py;
		pz = a.pz() + pz;
		energy = a.e() + energy;
		cachePhiRapidity();
	}

	/**
	 * Add to this particle another particle. Also add index of the added
	 * particle, which will be stored as an array. The method also precomputes
	 * rapidity, phi and et2 for fast retrival using getPhi, getRapidity
	 * methods.
	 * 
	 * @param a
	 * @param index
	 *            index of the particle to be stored
	 */
	public void add(ParticleD a, int index) {
		px = a.px() + px;
		py = a.py() + py;
		pz = a.pz() + pz;
		energy = a.e() + energy;
		cachePhiRapidity();
		consts.add(new Integer(index));
	}

	/**
	 * Get a hash code
	 */
	public int hashCode() {
		return hashCode() + (int) Double.doubleToRawLongBits(energy);
	}

}
