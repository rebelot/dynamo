use crate::topology::{atom::Atom, Topology};

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
    pub fn calc(&mut self, top: &mut Topology) {
        self.bonded.calc(top);
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
    pub fn calc(&mut self, top: &mut Topology) {
        for bond in self.bonds.iter_mut() {
            bond.calc(top);
        }
        for angle in self.angles.iter_mut() {
            angle.calc(top);
        }
        for dih in self.dihedrals.iter_mut() {
            dih.calc(top);
        }
    }
}

#[derive(Debug)]
pub struct Bond {
    k: f64,
    r0: f64,
    atoms: [usize; 2],
    disp_vec: [f64; DIM],
}

impl Bond {
    pub fn new(k: f64, r0: f64, atoms: [usize; 2]) -> Self {
        Self {
            k,
            r0,
            atoms,
            disp_vec: [0.0; DIM],
        }
    }

    fn calc(&mut self, top: &mut Topology) {
        let rvec = &mut self.disp_vec;

        let avec = &top.atoms[self.atoms[0]].pos;
        let bvec = &top.atoms[self.atoms[1]].pos;

        displace_vec(avec, bvec, rvec);
        let r2 = dot(rvec, rvec);
        let r = r2.sqrt();

        let mut u: f64 = 0.0;
        let mut f: f64 = 0.0;
        harmonic(&self.k, &self.r0, &r, &mut u, &mut f);
        f *= r2.sqrt().recip(); // shoudl optimize for fast invsqrt?

        let mut fbond: f64;
        for i in 0..DIM {
            fbond = f * self.disp_vec[i];
            top.atoms[self.atoms[0]].force[i] += fbond;
            top.atoms[self.atoms[1]].force[i] -= fbond;
        }
    }
}

#[derive(Debug)]
pub struct Angle {
    k: f64,
    t0: f64,
    atoms: [usize; 3],
    disp_vec: [f64; DIM],
}

impl Angle {
    pub fn new(k: f64, t0: f64, atoms: [usize; 3]) -> Self {
        Self {
            k,
            t0,
            atoms,
            disp_vec: [0.0; DIM],
        }
    }
    fn calc(&mut self, top: &mut Topology) {}

}

#[derive(Debug)]
pub struct Dihedral {
    k: f64,
    t0: f64,
    atoms: [usize; 4],
    disp_vec: [f64; DIM],
}

impl Dihedral {
    pub fn new(k: f64, t0: f64, atoms: [usize; 4]) -> Self {
        Self {
            k,
            t0,
            atoms,
            disp_vec: [0.0; DIM],
        }
    }
    fn calc(&mut self, top: &mut Topology) {}
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

fn dot(avec: &[f64; DIM], bvec: &[f64; DIM]) -> f64 {
    return avec.iter().zip(bvec.iter()).map(|(a, b)| a * b).sum();
}

/// Calculate the displacement vector between two points
fn displace_vec(avec: &[f64; DIM], bvec: &[f64; DIM], outvec: &mut [f64; DIM]) {
    for i in 0..DIM {
        outvec[i] = avec[i] - bvec[i];
    }
}

fn harmonic(k: &f64, r0: &f64, r: &f64, u: &mut f64, f: &mut f64) {
    let dr = r - r0;
    let kdr = k * dr;
    *u = 0.5 * kdr * dr;
    *f = -kdr;
}
