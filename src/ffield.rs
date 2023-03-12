use std::{cell::RefCell, rc::Rc};

use crate::topology::atom::Atom;

const DIM: usize = 3;

#[derive(Debug)]
pub struct FFList {
    pub nonbond: Nonbond,
    pub bond: Vec<Bond>,
    pub angle: Vec<Angle>,
    pub dihedral: Vec<Dihedral>,
}

impl FFList {
    pub fn new() -> FFList {
        FFList {
            nonbond: Nonbond::new(),
            bond: Vec::new(),
            angle: Vec::new(),
            dihedral: Vec::new(),
        }
    }
    pub fn calc(&mut self) {
        for bond in self.bond.iter_mut() {
            bond.calc();
        }
    }
}

#[derive(Debug)]
pub struct Bond {
    k: f64,
    r0: f64,
    atoms: [Rc<RefCell<Atom>>; 2],
    disp_vec: [f64; DIM],
}

impl Bond {
    pub fn new(k: f64, r0: f64, atoms: [Rc<RefCell<Atom>>; 2]) -> Self {
        Self {
            k,
            r0,
            atoms,
            disp_vec: [0.0; DIM],
        }
    }

    fn partial(&self, force: &f64) {
        let mut fbond: f64;
        for i in 0..DIM {
            fbond = force * self.disp_vec[i];
            self.atoms[0].borrow_mut().force[i] += fbond;
            self.atoms[1].borrow_mut().force[i] -= fbond;
        }
    }

    fn calc(&mut self) {
        let ai = self.atoms[0].borrow_mut();
        let aj = self.atoms[1].borrow_mut();
        let rvec = &mut self.disp_vec;

        displace_vec(&ai.pos, &aj.pos, rvec);
        let r2 = dot(rvec, rvec);
        let r = r2.sqrt();

        let mut u: f64 = 0.0;
        let mut f: f64 = 0.0;
        harmonic(&self.k, &self.r0, &r, &mut u, &mut f);
        f *= r2.sqrt().recip(); // shoudl optimize for fast invsqrt?
        self.partial(&f)
    }
}

#[derive(Debug)]
pub struct Angle {
    k: f64,
    t0: f64,
    atoms: [Rc<RefCell<Atom>>; 3],
}

impl Angle {
    pub fn new(k: f64, t0: f64, atoms: [Rc<RefCell<Atom>>; 3]) -> Self {
        Self { k, t0, atoms }
    }
}

#[derive(Debug)]
pub struct Dihedral {
    k: f64,
    t0: f64,
    atoms: [Rc<RefCell<Atom>>; 4],
}

impl Dihedral {
    pub fn new(k: f64, t0: f64, atoms: [Rc<RefCell<Atom>>; 4]) -> Self {
        Self { k, t0, atoms }
    }
}

#[derive(Debug)]
pub struct Nonbond {
    pub pairs: Vec<LJPair>,
}
impl  Nonbond {
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
    pub atoms: [Rc<RefCell<Atom>>; 2],
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

fn harmonic(k: &f64, r0: &f64, r: &f64, u: &mut f64, f: &mut f64 ) {
    let dr = r - r0;
    let kdr = k * dr;
    *u = 0.5 * kdr * dr;
    *f = -kdr;
}
