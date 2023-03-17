use crate::linalg::*;
use crate::topology::atom::Atom;

const DIM: usize = 3;

#[derive(Debug)]
pub struct Forces {
    pub bonded: BondedForces,
    pub nonbond: NonBondedPairs,
}
impl Forces {
    pub fn new() -> Forces {
        Forces {
            bonded: BondedForces::new(),
            nonbond: NonBondedPairs::new(),
        }
    }
    pub fn calc(&mut self, atoms: &mut Vec<Atom>) {
        self.bonded.calc(atoms);
    }
}

#[derive(Debug)]
pub struct BondedForces {
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub dihedrals: Vec<Dihedral>,
}

impl BondedForces {
    pub fn new() -> BondedForces {
        BondedForces {
            bonds: Vec::new(),
            angles: Vec::new(),
            dihedrals: Vec::new(),
        }
    }
    fn calc(&mut self, atoms: &mut Vec<Atom>) {
        for bond in self.bonds.iter_mut() {
            bond.calc(atoms);
        }
        for angle in self.angles.iter_mut() {
            angle.calc(atoms);
        }
        for dih in self.dihedrals.iter_mut() {
            dih.calc(atoms);
        }
    }
}

#[derive(Debug)]
pub struct Bond {
    k: f64,
    r0: f64,
    atoms: [usize; 2],
}

impl Bond {
    pub fn new(k: f64, r0: f64, atoms: [usize; 2]) -> Self {
        Self {
            k,
            r0,
            atoms,
        }
    }

    fn calc(&mut self, atoms: &mut Vec<Atom>) {

        let ri = &atoms[self.atoms[0]].pos;
        let rj = &atoms[self.atoms[1]].pos;

        let rij = displace_vec(ri, rj);
        let norm2_rij = norm2(&rij);
        let norm_rij = norm2_rij.sqrt();

        let u = &mut 0.0;
        let f = &mut 0.0;
        harmonic(&self.k, &self.r0, &norm_rij, u, f);
        *f *= norm_rij.recip();
        // *f *= norm2_rij.sqrt().recip(); // will this be faster?

        let di = &mut 0.0;
        for i in 0..DIM {
            *di = *f * rij[i];
            atoms[self.atoms[0]].force[i] += *di;
            atoms[self.atoms[1]].force[i] -= *di;
        }
    }
}

#[derive(Debug)]
pub struct Angle {
    k: f64,
    t0: f64,
    atoms: [usize; 3],
}

impl Angle {
    pub fn new(k: f64, t0: f64, atoms: [usize; 3]) -> Self {
        Self {
            k,
            t0,
            atoms,
        }
    }
    fn calc(&mut self, atoms: &mut Vec<Atom>) {
        let ri = &atoms[self.atoms[0]].pos;
        let rj = &atoms[self.atoms[1]].pos;
        let rk = &atoms[self.atoms[2]].pos;

        let rij = displace_vec(ri, rj);
        let rkj = displace_vec(rk, rj);

        let norm_rij_1 = norm2(&rij).sqrt().recip();
        let norm_rkj_1 = norm2(&rkj).sqrt().recip();

        let norm_rij_rkj_1 = norm_rij_1 * norm_rkj_1;

        let cos_t = dot(&rij, &rkj) * norm_rij_rkj_1;
        let t = cos_t.acos();

        let u = &mut 0.0;
        let f = &mut 0.0;

        harmonic(&self.k, &self.t0, &t, u, f);

        *f *= -(1.0 - cos_t * cos_t).sqrt().recip() * norm_rij_rkj_1;
        let (di, dj, dk) = (&mut 0.0, &mut 0.0, &mut 0.0);
        for i in 0..DIM {
            *di = -rkj[i];
            *dk = -rij[i];
            *dj = -*di - *dk;
            atoms[self.atoms[0]].force[i] += *f * *di;
            atoms[self.atoms[1]].force[i] += *f * *dj;
            atoms[self.atoms[2]].force[i] += *f * *dk;
        }
    }
}

#[derive(Debug)]
pub struct Dihedral {
    k: f64,
    t0: f64,
    atoms: [usize; 4],
}

impl Dihedral {
    pub fn new(k: f64, t0: f64, atoms: [usize; 4]) -> Self {
        Self {
            k,
            t0,
            atoms,
        }
    }
    fn calc(&mut self, atoms: &mut Vec<Atom>) {
        let ri = &atoms[self.atoms[0]].pos;
        let rj = &atoms[self.atoms[1]].pos;
        let rk = &atoms[self.atoms[2]].pos;
        let rl = &atoms[self.atoms[3]].pos;

        let rij = displace_vec(ri, rj);
        let rjk = displace_vec(rj, rk);
        let rkl = displace_vec(rk, rl);

        let nijk = cross(&rij, &rjk);
        let njkl = cross(&rjk, &rkl);

        let t = angle(&nijk, &njkl);

        let u = &mut 0.0;
        let f = &mut 0.0;

        harmonic(&self.k, &self.t0, &t, u, f);

    }
}

#[derive(Debug)]
pub struct NonBondedPairs {
    pub pairs: Vec<LJPair>,
}

impl NonBondedPairs {
    pub fn new() -> Self {
        Self { pairs: Vec::new() }
    }
}

#[derive(Debug)]
pub struct LJParams {
    pub c12: f64,
    pub c6: f64,
}
impl LJParams {
    pub fn new(c12: f64, c6: f64) -> Self {
        Self { c12, c6 }
    }
}

#[derive(Debug)]
pub struct LJPair {
    pub atoms: [usize; 2],
}

fn harmonic(k: &f64, r0: &f64, r: &f64, u: &mut f64, f: &mut f64) {
    let dr = r - r0;
    let kdr = k * dr;
    *u = 0.5 * kdr * dr;
    *f = -kdr;
}
