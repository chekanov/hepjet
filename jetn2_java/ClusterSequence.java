import java.util.*;
import java.io.*;

/**
 * Main class for jet clustering with generalized kt algorithms.
 * Implements optimizations relevant to HEP events:
 *   1. Geometric factorization -- reduces overall complexity to O(N^2).
 *      http://arxiv.org/pdf/hep-ph/0512210.pdf
 *      http://www.lpthe.jussieu.fr/~salam/repository/docs/kt-cgta-v2.pdf
 *   2. Tiling -- optimizes the number of pairwise distances to compute.
 *   3. Linked list -- optimizes insertion and deletion, O(1) complexity.
 * Contains references to the first PseudeJet and y-phi TileGrid.
 * Safe to only use one ClusterSequence instance per thread.
 *
 * @author Ivan Pogrebnyak
 *
 */
class ClusterSequence {

  private final static double twopi = 2*Math.PI;
  private final static double max_rap = 1e5; // From FastJet
  private static double sq(double x) { return x*x; }

  private byte   alg;
  private double jetR2;
  private int    num;

  private PseudoJet first;
  private TileGrid  grid;
  private boolean   use_grid;

  /**
   * ClusterSequence constructor.
   * @param alg Select jet clustering algorithm. Accepted values are:
   *         "kt": $d_{ij} = \min(k_{ti}^{-2},k_{tj}^{-2})\frac{\Delta^2_{ij}}{R^2}$
   *     "antikt": $d_{ij} = \min(k_{ti}^{2},k_{tj}^{2})\frac{\Delta^2_{ij}}{R^2}$
   *  "cambridge" or "ca": $d_{ij} = \frac{\Delta^2_{ij}}{R^2}$
   * @param jetR Algorithm radius parameter, R.
   */
  public ClusterSequence(String alg, double jetR) {
    if (alg.equalsIgnoreCase("kt")) this.alg = 1;
    else if (alg.equalsIgnoreCase("antikt")) this.alg = -1;
    else if (alg.equalsIgnoreCase("ca")) this.alg = 0;
    else if (alg.equalsIgnoreCase("cambridge")) this.alg = 0;
    else {
      System.out.println(
        "Warning: Unrecognized clustering algorithm: "+alg+
        ".\nDefaulting to antikt."
      );
      this.alg = -1;
    }
    this.jetR2 = jetR*jetR;

    grid = new TileGrid(jetR,5);
  }

  /**
   * private PseudoJet class.
   * Combined 4-momentum and linked list node functionalities.
   * Represents intermediate clustering sequence pseudo-jets.
   */
  private class PseudoJet {
    public double px, py, pz, E, rap, phi, Rij, diB, dij;
    public int id;
    public PseudoJet prev, next, near;
    public Tile t;
    public PseudoJet tprev, tnext;
    public ArrayList<Integer> consts;

    public PseudoJet(double px, double py, double pz, double E) {
      this.px = px;
      this.py = py;
      this.pz = pz;
      this.E  = E;
      this.id = num++;

      final double pt2 = px*px + py*py;
      final double abs_pz = (pz < 0. ? -pz : pz);

      phi = (pt2 == 0. ? 0. : Math.atan2(py,px)) + Math.PI;
      if (phi >= twopi) phi -= twopi;
      else if (phi < 0.) phi += twopi;

      // Replicated FastJet rapidity calculation
      // for compatibility in cases of unphysical 4-momenta
      if (E == abs_pz && pt2 == 0.) {
        // Point has infinite rapidity -- convert that into a very large
        // number, but in such a way that different 0-pt momenta will have
        // different rapidities (so as to lift the degeneracy between
        // them) [this can be relevant at parton-level]
        rap = max_rap + abs_pz;
        if (pz < 0.) rap = -rap;
      } else {
        // get the rapidity in a way that's modestly insensitive to roundoff
        // error when things pz,E are large (actually the best we can do without
        // explicit knowledge of mass) and force non tachyonic mass
        double m2_pt2 = (E+pz)*(E-pz);
        if (m2_pt2 < pt2) m2_pt2 = pt2;
        rap = 0.5*Math.log(m2_pt2/sq(E+abs_pz));
        if (pz > 0.) rap = -rap;
      }

      switch (alg) {
        case -1: if (pt2==0.) diB = Double.MAX_VALUE;
                 else diB = 1./pt2; // antikt
                 break;
        case  1: diB = pt2; break;  // kt
        case  0: diB = 1.;  break;  // cambridge
      }

      Rij = Double.MAX_VALUE;
    }

    public void remove() {
      if (prev==null) first = next;
      else prev.next = next;
      if (next!=null) next.prev = prev;

      if (use_grid) {
        if (tprev==null) t.first = tnext;
        else tprev.tnext = tnext;
        if (tnext!=null) tnext.tprev = tprev;
      }
    }

    public void merge() {
      // add 4-momenta
      first.prev = new PseudoJet(
        px + near.px, py + near.py,
        pz + near.pz, E  + near.E );
      first.prev.next = first;
      first = first.prev;

      // merge constituents
      if (this.consts!=null || near.consts!=null) {
        if (this.consts!=null) {
          first.consts = this.consts;
          if (near.consts!=null) first.consts.addAll(near.consts);
          else first.consts.add(near.id);
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

    public boolean update_near(PseudoJet p, boolean both) {
      double deltaPhi = Math.abs(phi-p.phi);
      if (deltaPhi > Math.PI) deltaPhi = twopi - deltaPhi;
      double Ril = sq(rap-p.rap) + sq(deltaPhi);

      if (Ril < Rij) { Rij = Ril; near = p; }
      if (both) {
        if (Ril < p.Rij) { p.Rij = Ril; p.near = this; return true; }
        else return false;
      } else return false;
    }

    public void update_dij() {
      if (near==null) dij = Double.MAX_VALUE;
      else dij = Math.min(diB,near.diB)*Rij/jetR2;
    }

    public void rm_near() {
      near = null;
      Rij  = Double.MAX_VALUE;
    }
  }

  /**
   * private Tile class. Implements grid tile.
   */
  class Tile {
    int irap, iphi;
    double rapc, phic; // center coordinates
    PseudoJet first;
    Tile(int irap, int iphi, double rapc, double phic) {
      this.irap = irap;
      this.iphi = iphi;
      this.rapc = rapc;
      this.phic = phic;
    }
  }

  /**
   * private TileGrid class
   */
  private class TileGrid {
    private final Tile[][] tiles;
    private final int nrap, nphi;
    private final double d, r;
    private final double max_rap;

    public TileGrid(double R, double max_rap) {
      nphi = (int)(twopi/R);
      d = twopi/nphi;
      r = d/2;
      nrap = 2*((int)(max_rap/r) + 1);
      this.max_rap = nrap*r;

      tiles = new Tile[nphi][nrap];
      for (int irap=0; irap<nrap; ++irap)
        for (int iphi=0; iphi<nphi; ++iphi)
          tiles[iphi][irap] = new Tile(
            irap, iphi,
            (irap-nrap/2)*d + r, iphi*d + r
          );
    }

    public void clear() {
      for (int irap=0; irap<nrap; ++irap)
        for (int iphi=0; iphi<nphi; ++iphi)
          tiles[iphi][irap].first = null;
    }

    public void add(PseudoJet p) {
      int irap = (int)((p.rap+max_rap)/d);
      if (irap < 0) irap = 0;
      else if (irap >= nrap) irap = nrap-1;

      p.t = tiles[ (int)(p.phi/d) ][ irap ];

      if (p.t.first==null) {
        p.t.first = p;
      } else {
        p.tnext = p.t.first;
        p.t.first.tprev = p;
        p.t.first = p;
      }
    }

    private void within_tile(PseudoJet p, Tile t, boolean both) {
      for (PseudoJet q=t.first; q!=null; q=q.tnext)
        if (p.update_near(q, both)) q.update_dij();
    }

    public void update_near(PseudoJet p, boolean both) {
      // own tile
      for (PseudoJet q=p.t.first; q!=null; q=q.tnext)
        if (q!=p)
          if (p.update_near(q, both)) q.update_dij();

      final boolean // lazy
        td = ( sq(p.phi-p.t.phic+r) < p.Rij ),
        tu = ( sq(p.phi-p.t.phic-r) < p.Rij ),
        tl = ( sq(p.rap-p.t.rapc+r) < p.Rij ),
        tr = ( sq(p.rap-p.t.rapc-r) < p.Rij );

      if (p.t.irap!=0) {
        if (tl) { within_tile(p, tiles[p.t.iphi][p.t.irap-1], both);
          if (td) within_tile(p, tiles[p.t.iphi==0 ? nphi-1 : p.t.iphi-1][p.t.irap-1], both);
          if (tu) within_tile(p, tiles[nphi-p.t.iphi==1 ? 0 : p.t.iphi+1][p.t.irap-1], both);
        }
      }

      if (td) within_tile(p, tiles[p.t.iphi==0 ? nphi-1 : p.t.iphi-1][p.t.irap], both);
      if (tu) within_tile(p, tiles[nphi-p.t.iphi==1 ? 0 : p.t.iphi+1][p.t.irap], both);

      if (nrap-p.t.irap!=1) {
        if (tr) { within_tile(p, tiles[p.t.iphi][p.t.irap+1], both);
          if (td) within_tile(p, tiles[p.t.iphi==0 ? nphi-1 : p.t.iphi-1][p.t.irap+1], both);
          if (tu) within_tile(p, tiles[nphi-p.t.iphi==1 ? 0 : p.t.iphi+1][p.t.irap+1], both);
        }
      }

    }
  }

  /**
   * Form jets from particles in a single event.
   * @param particles List of input particles.
   * @param min_jet_pt Save output jets only with pT >= min_jet_pt.
   * @return Clustered jets.
   */
  public ArrayList<ParticleD> cluster(List<ParticleD> particles, double min_jet_pt) {
    int n = particles.size();
    num = 0; // start assigning PseudoJet id from 0

    // initialize the grid
    use_grid = (n>50);

    ArrayList<ParticleD> jets = new ArrayList<ParticleD>();
    PseudoJet p;

    if (n==0) return jets;

    // read input particles -------------------------------
    first = new PseudoJet(
      particles.get(0).px(), particles.get(0).py(),
      particles.get(0).pz(), particles.get(0).e ()
    );
    p = first;
    if (use_grid) grid.add(p);
    for (int i=1; i<n; ++i) {
      p.next = new PseudoJet(
        particles.get(i).px(), particles.get(i).py(),
        particles.get(i).pz(), particles.get(i).e ()
      );
      p.next.prev = p;
      p = p.next;
      if (use_grid) grid.add(p);
    }

    // find original nearest neighbors --------------------
    if (use_grid) { // using grid

      for (p=first; p!=null; p=p.next)
        grid.update_near(p,false);

    } else { // no grid

      for (p=first.next; p!=null; p=p.next)
        for (PseudoJet q=first; q!=p; q=q.next)
          p.update_near(q,true);

    }

    // calculate minimum pairwise distances ---------------
    for (p=first; p!=null; p=p.next) p.update_dij();

    // loop until PseudoJets are used up ------------------
    while (first != null) {

      if (n<50) {
        use_grid = false;
        grid.clear();
      }

      p = first;
      double dist = p.diB;
      boolean merge = false;

      // find smallest distance
      for (PseudoJet q=first; q!=null; q=q.next) {
        if (q.dij < dist) { p = q; dist = q.dij; merge = true; }
      }

      if (p.Rij > jetR2) {
        for (PseudoJet q=first.next; q!=null; q=q.next) {
          if (q.diB < dist) { p = q; dist = q.diB; merge = false; }
        }
      }

      // Either merge or identify a jet
      if (merge) {

        // merge particles
        p.merge();

        // the new particle is first
        if (use_grid) grid.add(first);

        // print clustering step
        // System.out.format("%4d & %4d | d = %.5e\n", p.id, p.near.id, dist);

        // recompute pairwise distances
        if (use_grid) { // using grid

          // for the new PseudoJet
          grid.update_near(first,true);
          first.update_dij();

          // for the rest
          for (PseudoJet p1=first.next; p1!=null; p1=p1.next) {
            if (p1.near==p || p1.near==p.near) {
              p1.rm_near();
              grid.update_near(p1,false);
              p1.update_dij();
            }
          }

        } else { // no grid

          // for the new PseudoJet
          for (PseudoJet q=first.next; q!=null; q=q.next)
            if ( first.update_near(q,true) ) q.update_dij();
          first.update_dij();

          // for the rest
          for (PseudoJet p1=first.next; p1!=null; p1=p1.next) {
            if (p1.near==p || p1.near==p.near) {
              p1.rm_near();
              for (PseudoJet p2=first; p2!=null; p2=p2.next) {
                if (p1!=p2) p1.update_near(p2,false);
              }
              p1.update_dij();
            }
          }

        }

      } else {
        // identify as jet
        if ( Math.sqrt(sq(p.px)+sq(p.py)) >= min_jet_pt ) {
          ParticleD jet = new ParticleD(p.px, p.py, p.pz, p.E);
          jets.add(jet);
          if (p.consts==null) jet.addConstituent(p.id);
          else jet.setConstituents(p.consts);
        }

        // print clustering step
        // System.out.format("%4d Jet    | d = %.5e\n", p.id, dist);

        p.remove(); // "remove"

        // recompute pairwise distances
        if (use_grid) { // using grid

          for (PseudoJet p1=first; p1!=null; p1=p1.next) {
            if (p1.near==p) {
              p1.rm_near();
              grid.update_near(p1,false);
              p1.update_dij();
            }
          }

        } else { // no grid

          for (PseudoJet p1=first; p1!=null; p1=p1.next) {
            if (p1.near==p) {
              p1.rm_near();
              for (PseudoJet p2=first; p2!=null; p2=p2.next) {
                if (p1!=p2) p1.update_near(p2,false);
              }
              p1.update_dij();
            }
          }

        }

      }

      --n;

    }

    // remove all particles --------------------------
    first = null;
    if (use_grid) grid.clear();

    return jets;
  }
}
