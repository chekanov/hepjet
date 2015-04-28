import java.util.*;
import java.io.*;

/**
 * Main class to perform kt and anti-kT jet algorithm using a geometrical
 * approach as implemented in the FastJet C++ algorithm.
 *
 * @author Ivan Pogrebnyak, Sergei Chekanov
 *
 */

class ClusterSequence {

	private final static double twopi = 2 * Math.PI;

	private static double sq(double x) {
		return x * x;
	}

	private byte alg;
	private double jetR2;
	private int num;

	private pseudoJet first;
	private tileGrid grid;
	private boolean use_grid;

	// ****************************************************************
	// Constructor ****************************************************
	public ClusterSequence(String alg, double jetR) {
		if (alg.equalsIgnoreCase("kt"))
			this.alg = 1;
		else if (alg.equalsIgnoreCase("antikt"))
			this.alg = 2;
		else if (alg.equalsIgnoreCase("ca"))
			this.alg = 3;
		else {
			System.out.println("Warning: Unrecognized clustering algorithm: "
					+ alg + ".\nDefaulting to antikt.");
			this.alg = 2;
		}
		this.jetR2 = jetR * jetR;

		grid = new tileGrid(jetR, 5);
	}

	// ****************************************************************
	// pseudoJet class ************************************************
	private class pseudoJet {
		public double px, py, pz, E, rap, phi, Rij, diB, dij;
		public int id;
		public pseudoJet prev, next, near;
		public tile t;
		public pseudoJet tprev, tnext;
		public ArrayList<Integer> consts;

		public pseudoJet(double px, double py, double pz, double E) {
			this.px = px;
			this.py = py;
			this.pz = pz;
			this.E = E;
			this.id = num++;

			if (E == pz)
				rap = Double.MAX_VALUE;
			else if (E == -pz)
				rap = -Double.MAX_VALUE;
			else
				rap = 0.5 * Math.log((E + pz) / (E - pz));
			phi = (px == 0. && py == 0. ? 0. : Math.atan2(py, px)) + Math.PI;

			switch (alg) {
			case 1:
				diB = pt2();
				break; // kt
			case 2:
				diB = 1. / pt2();
				if (Double.isInfinite(diB)) diB = Double.MAX_VALUE;
				break; // antikt
			case 3:
				diB = 1.;
				break; // cambridge
			}
			Rij = Double.MAX_VALUE;
		}

		public void remove() {
			if (prev == null)
				first = next;
			else
				prev.next = next;
			if (next != null)
				next.prev = prev;

			if (use_grid) {
				if (tprev == null)
					t.first = tnext;
				else
					tprev.tnext = tnext;
				if (tnext != null)
					tnext.tprev = tprev;
			}
		}

		public void merge() {
			// add 4-momenta
			first.prev = new pseudoJet(px + near.px, py + near.py,
					pz + near.pz, E + near.E);
			first.prev.next = first;
			first = first.prev;

			// merge constituents
			if (this.consts != null || near.consts != null) {
				if (this.consts != null) {
					first.consts = this.consts;
					if (near.consts != null)
						first.consts.addAll(near.consts);
					else
						first.consts.add(near.id);
				} else {
					first.consts = near.consts;
					first.consts.add(this.id);
				}
			} else {
				first.consts = new ArrayList<Integer>();
				first.consts.add(this.id);
				first.consts.add(near.id);
			}

			// remove parent particles
			this.remove();
			near.remove();
		}

		public double pt2() {
			return px * px + py * py;
		}

		public boolean update_near(pseudoJet p, boolean both) {
			double deltaPhi = Math.abs(phi - p.phi);
			if (deltaPhi > Math.PI)
				deltaPhi = twopi - deltaPhi;
			double Ril = sq(rap - p.rap) + sq(deltaPhi);

			if (Ril < Rij) {
				Rij = Ril;
				near = p;
			}
			if (both) {
				if (Ril < p.Rij) {
					p.Rij = Ril;
					p.near = this;
					return true;
				} else
					return false;
			} else
				return false;
		}

		public void update_dij() {
			if (near == null)
				dij = Double.MAX_VALUE;
			else
				dij = Math.min(diB, near.diB) * Rij / jetR2;
		}

		public void rm_near() {
			near = null;
			Rij = Double.MAX_VALUE;
		}
	}

	// ****************************************************************
	// tile class *****************************************************
	private class tile {
		int irap, iphi;
		double rapc, phic; // center coordinates
		pseudoJet first;

		tile(int irap, int iphi, double rapc, double phic) {
			this.irap = irap;
			this.iphi = iphi;
			this.rapc = rapc;
			this.phic = phic;
		}
	}

	// ****************************************************************
	// tileGrid class *************************************************
	private class tileGrid {
		private final tile[][] tiles;
		private final int nrap, nphi;
		private final double d, r;
		private final double max_rap;

		public tileGrid(double R, double max_rap) {
			nphi = (int) (twopi / R);
			d = twopi / nphi;
			r = d / 2;
			nrap = 2 * ((int) (max_rap / r) + 1);
			this.max_rap = nrap * r;

			tiles = new tile[nphi][nrap];
			for (int irap = 0; irap < nrap; ++irap)
				for (int iphi = 0; iphi < nphi; ++iphi)
					tiles[iphi][irap] = new tile(irap, iphi, (irap - nrap / 2)
							* d + r, iphi * d + r);
		}

		public void clear() {
			for (int irap = 0; irap < nrap; ++irap)
				for (int iphi = 0; iphi < nphi; ++iphi)
					tiles[iphi][irap].first = null;
		}

		public void add(pseudoJet p) {
			int irap = (int) ((p.rap + max_rap) / d);
			if (irap < 0)
				irap = 0;
			else if (irap >= nrap)
				irap = nrap - 1;

			p.t = tiles[(int) (p.phi / d)][irap];

			if (p.t.first == null) {
				p.t.first = p;
			} else {
				p.tnext = p.t.first;
				p.t.first.tprev = p;
				p.t.first = p;
			}
		}

		private void within_tile(pseudoJet p, tile t, boolean both) {
			for (pseudoJet q = t.first; q != null; q = q.tnext)
				if (p.update_near(q, both))
					q.update_dij();
		}

		public void update_near(pseudoJet p, boolean both) {
			// own tile
			for (pseudoJet q = p.t.first; q != null; q = q.tnext)
				if (q != p)
					if (p.update_near(q, both))
						q.update_dij();

			final boolean // lazy
			tl = (sq(p.phi - p.t.phic + r) < p.Rij), tr = (sq(p.phi - p.t.phic
					- r) < p.Rij), td = (sq(p.rap - p.t.rapc + r) < p.Rij), tu = (sq(p.rap
					- p.t.rapc - r) < p.Rij);

			if (p.t.irap != 0) {
				if (tl) {
					within_tile(p, tiles[p.t.iphi][p.t.irap - 1], both);
					if (tu)
						within_tile(p, tiles[p.t.iphi == 0 ? nphi - 1
								: p.t.iphi - 1][p.t.irap - 1], both);
					if (td)
						within_tile(p, tiles[nphi - p.t.iphi == 1 ? 0
								: p.t.iphi + 1][p.t.irap - 1], both);
				}
			}

			if (tu)
				within_tile(
						p,
						tiles[p.t.iphi == 0 ? nphi - 1 : p.t.iphi - 1][p.t.irap],
						both);
			if (td)
				within_tile(
						p,
						tiles[nphi - p.t.iphi == 1 ? 0 : p.t.iphi + 1][p.t.irap],
						both);

			if (nrap - p.t.irap != 1) {
				if (tr) {
					within_tile(p, tiles[p.t.iphi][p.t.irap + 1], both);
					if (tu)
						within_tile(p, tiles[p.t.iphi == 0 ? nphi - 1
								: p.t.iphi - 1][p.t.irap + 1], both);
					if (td)
						within_tile(p, tiles[nphi - p.t.iphi == 1 ? 0
								: p.t.iphi + 1][p.t.irap + 1], both);
				}
			}

		}
	}


	/**
	 * Main class to perform a clustering
	 * @param particles list with input particles
	 * @param minpt min pT for jets
	 * @return output jets above certain pT
	 */
	public ArrayList<ParticleD> cluster(List<ParticleD> particles, double minpt) {
		final int n = particles.size();
		num = 0; // start assigning pseudoJet id from 0
			
		// initialize the grid
		use_grid = (n > 50);

		ArrayList<ParticleD> jets = new ArrayList<ParticleD>();
		pseudoJet p;

		if (n == 0)
			return jets;

		// read input particles -------------------------------
		p = first = new pseudoJet(particles.get(0).px(), particles.get(0).py(),
				particles.get(0).pz(), particles.get(0).e());
		if (use_grid)
			grid.add(p);
		for (int i = 1; i < n; ++i) {
			p.next = new pseudoJet(particles.get(i).px(),
					particles.get(i).py(), particles.get(i).pz(), particles
							.get(i).e());
			p.next.prev = p;
			p = p.next;
			if (use_grid)
				grid.add(p);
		}

		// find original nearest neighbors --------------------
		if (use_grid) { // using grid

			for (p = first; p != null; p = p.next)
				grid.update_near(p, false);

		} else { // no grid

			for (p = first.next; p != null; p = p.next)
				for (pseudoJet q = first; q != p; q = q.next)
					p.update_near(q, true);

		}

		// calculate minimum pairwise distances ---------------
		for (p = first; p != null; p = p.next)
			p.update_dij();

		// loop until pseudoJets are used up ------------------
		while (first != null) {

			p = first;
			double dist = p.diB;
			boolean merge = false;

			// find smallest distance
			for (pseudoJet q = first; q != null; q = q.next) {
				if (q.dij < dist) {
					p = q;
					dist = q.dij;
					merge = true;
				}
			}

			if (p.Rij > jetR2) {
				for (pseudoJet q = first.next; q != null; q = q.next) {
					if (q.diB < dist) {
						p = q;
						dist = q.diB;
						merge = false;
					}
				}
			}

			// Either merge or identify a jet
			if (merge) {

				// merge particles
				p.merge();

				// the new particle is first
				if (use_grid)
					grid.add(first);

				// print clustering step
				// System.out.format("%4d & %4d | d = %.5e\n", p.id, p.near.id, dist);

				// recompute pairwise distances
				if (use_grid) { // using grid

					// for the new pseudoJet
					grid.update_near(first, true);
					first.update_dij();

					// for the rest
					for (pseudoJet p1 = first.next; p1 != null; p1 = p1.next) {
						if (p1.near == p || p1.near == p.near) {
							p1.rm_near();
							grid.update_near(p1, false);
							p1.update_dij();
						}
					}

				} else { // no grid

					// for the new pseudoJet
					for (pseudoJet q = first.next; q != null; q = q.next)
						if (first.update_near(q, true))
							q.update_dij();
					first.update_dij();

					// for the rest
					for (pseudoJet p1 = first.next; p1 != null; p1 = p1.next) {
						if (p1.near == p || p1.near == p.near) {
							p1.rm_near();
							for (pseudoJet p2 = first; p2 != null; p2 = p2.next) {
								if (p1 != p2)
									p1.update_near(p2, false);
							}
							p1.update_dij();
						}
					}

				}

			} else {
				// identify as jet
				if (Math.sqrt(p.pt2()) > minpt) {
					ParticleD jee = new ParticleD(p.px, p.py, p.pz, p.E);
					jets.add(jee);
					if (p.consts == null)
						jee.addConstituent(p.id);
					else
						jee.setConstituents(p.consts);
				}

				// print clustering step
				// System.out.format("%4d Jet    | d = %.5e\n", p.id, dist);

				// "remove"
				p.remove();

				// recompute pairwise distances
				if (use_grid) { // using grid

					for (pseudoJet p1 = first; p1 != null; p1 = p1.next) {
						if (p1.near == p) {
							p1.rm_near();
							grid.update_near(p1, false);
							p1.update_dij();
						}
					}

				} else { // no grid

					for (pseudoJet p1 = first; p1 != null; p1 = p1.next) {
						if (p1.near == p) {
							p1.rm_near();
							for (pseudoJet p2 = first; p2 != null; p2 = p2.next) {
								if (p1 != p2)
									p1.update_near(p2, false);
							}
							p1.update_dij();
						}
					}

				}

			}

		}

		// remove all particles --------------------------
		first = null;
		if (use_grid)
			grid.clear();

		return jets;
	}
}
