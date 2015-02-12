// A simple particle class to keep a particle or jet
// S.Chekanov

#include "ParticleD.h"

ParticleD::ParticleD()
{
	consts = std::vector<int>();
        consts.push_back(-1); // itself
}

ParticleD::ParticleD(double m_px, double m_py, double m_pz, double energy)
{
	this->m_px = m_px;
	this->m_py = m_py;
	this->m_pz = m_pz;
	this->energy = energy;
	cachePhiRapidity();
	consts = std::vector<int>();
        consts.push_back(-1); // itself
}

void ParticleD::setPxPyPzE(double m_px, double m_py, double m_pz, double energy)
{
	this->m_px = m_px;
	this->m_py = m_py;
	this->m_pz = m_pz;
	this->energy = energy;
	cachePhiRapidity();
}

double ParticleD::eta()
{
	if (m_px == 0.0 && m_py == 0.0)
	{
		return -999;
	}
	if (m_pz == 0.0)
	{
		return -999;
	}
	double pt2 = m_px*m_px + m_py*m_py;
	double theta = atan2(sqrt(pt2),m_pz);
	if (theta < 0)
	{
		theta += PI2;
	}
	return -1*log(tan(theta / 2));
}

double ParticleD::et2()
{
      double pt2x = perp2();
      double m_et2 = pt2x == 0 ? 0 : e()*e() * pt2x / (pt2x + pz()*pz());
      return m_et2;
}

double ParticleD::et()
{
	double etet = et2();
	return e() < 0.0 ? - sqrt(etet) : sqrt(etet);
}

double ParticleD::rapidity()
{
	m_rapidity = -10e10;
	if (energy > pz())
	{
		m_rapidity = 0.5*log((energy + m_pz) / (energy - m_pz));
	}
	return m_rapidity;
}

double ParticleD::mag()
{
	return sqrt(m_px*m_px + m_py*m_py + m_pz*m_pz);
}

double ParticleD::perp2()
{
        m_pt2=m_px*m_px + m_py*m_py;
	return m_pt2;
}

double ParticleD::perp()
{
	return sqrt(perp2());
}

void ParticleD::setEnergy(double energy)
{
	this->energy = energy;

}

int ParticleD::compareTo(ParticleD *o)
{
	if (getPt2() < o->getPt2())
	{
		return 1;
	}
	if (getPt2() > o->getPt2())
	{
		return -1;
	}
	return 0;
}

ParticleD *ParticleD::copy(ParticleD *o)
{
	ParticleD *tmp = new ParticleD(m_px,m_py,m_pz,energy);
	tmp->setRapidity(m_rapidity);
	tmp->setPhi(m_phi);
	tmp->setPt2(m_pt2);
	return tmp;
}


double ParticleD::px()
{
	return m_px;
}

double ParticleD::py()
{
	return m_py;
}

double ParticleD::pz()
{
	return m_pz;
}

double ParticleD::getRapidity()
{
	return m_rapidity;
}

double ParticleD::getPt2()
{
	return m_pt2;
}

double ParticleD::getPt()
{
	return sqrt(m_pt2);
}

double ParticleD::getPhi()
{
	return m_phi;
}

std::vector<int> ParticleD::getConstituents()
{
	return consts;
}


double ParticleD::phi()
{
	if (m_px == 0)
	{
		return 0.0;
	}
	m_phi = atan2(m_py,m_px);
	//if (m_phi < 0)
        //  	m_phi = PI2 + m_phi;
	return m_phi;
}

double ParticleD::e()
{
	return energy;
}

void ParticleD::setPx(double m_px)
{
	this->m_px = m_px;
}

void ParticleD::setPy(double m_py)
{
	this->m_py = m_py;
}

void ParticleD::setPz(double m_pz)
{
	this->m_pz = m_pz;
}

void ParticleD::setRapidity(double m_rapidity)
{
	this->m_rapidity = m_rapidity;
}

void ParticleD::setPhi(double m_phi)
{
	this->m_phi = m_phi;
}

void ParticleD::setPt2(double m_pt2)
{
	this->m_pt2 = m_pt2;
}

void ParticleD::cachePhiRapidity()
{
	m_rapidity = rapidity();
	m_phi = phi();
	m_pt2 = perp2();
}

void ParticleD::add(ParticleD *a)
{
	m_px = a->px() + m_px;
	m_py = a->py() + m_py;
	m_pz = a->pz() + m_pz;
	energy = a->e() + energy;
	cachePhiRapidity();
}

void ParticleD::add(ParticleD *a, int index)
{
	m_px = a->px() + m_px;
	m_py = a->py() + m_py;
	m_pz = a->pz() + m_pz;
	energy = a->e() + energy;
	cachePhiRapidity();
	consts.push_back(int(index));
}

bool ParticleD::operator<(ParticleD& score) const
{
  if(m_pt2 < score.getPt2())
    return true;      
  else if (m_pt2 == score.getPt2())
    return true;
  else return false;
}

