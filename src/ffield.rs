use crate::topology::atom::Atom;
use crate::{Rvec, DIM};

mod functions;

pub struct Forces {
    pub bonds: Vec<Box<dyn BondedInteraction>>,
    pub angles: Vec<Box<dyn AngleInteraction>>,
    pub torsions: Vec<Box<dyn TorsionInteraction>>,
}

impl Forces {
    pub fn new() -> Forces {
        Forces {
            bonds: Vec::new(),
            angles: Vec::new(),
            torsions: Vec::new(),
        }
    }

    pub fn calc(&mut self, atoms: &mut Vec<Atom>) {
        let mut u_tot = 0.0;
        self.bonds.iter().for_each(|bond| {
            let [a1, a2] = bond.atoms();
            let (u, f) = bond.calc(&atoms[a1].pos, &atoms[a2].pos);

            for i in 0..DIM {
                atoms[a1].force[i] += f[0][i];
                atoms[a2].force[i] += f[1][i];
            }
            u_tot += u;
        });
        self.angles.iter().for_each(|angle| {
            let [a1, a2, a3] = angle.atoms();
            let (u, f) = angle.calc(&atoms[a1].pos, &atoms[a2].pos, &atoms[a3].pos);
            for i in 0..DIM {
                atoms[a1].force[i] += f[0][i];
                atoms[a2].force[i] += f[1][i];
                atoms[a3].force[i] += f[2][i];
            }
            u_tot += u;
        });
        self.torsions.iter().for_each(|torsion| {
            let [a1, a2, a3, a4] = torsion.atoms();
            let (u, f) = torsion.calc(
                &atoms[a1].pos,
                &atoms[a2].pos,
                &atoms[a3].pos,
                &atoms[a4].pos,
            );
            for i in 0..DIM {
                atoms[a1].force[i] += f[0][i];
                atoms[a2].force[i] += f[1][i];
                atoms[a3].force[i] += f[2][i];
                atoms[a4].force[i] += f[3][i];
            }
            u_tot += u;
        });
    }
}

pub trait BondedInteraction {
    fn calc(&self, ri: &Rvec, rj: &Rvec) -> (f32, [Rvec; 2]);
    fn atoms(&self) -> [usize; 2];
}

pub trait AngleInteraction {
    fn calc(&self, ri: &Rvec, rj: &Rvec, rk: &Rvec) -> (f32, [Rvec; 3]);
    fn atoms(&self) -> [usize; 3];
}

pub trait TorsionInteraction {
    fn calc(&self, ri: &Rvec, rj: &Rvec, rk: &Rvec, rl: &Rvec) -> (f32, [Rvec; 4]);
    fn atoms(&self) -> [usize; 4];
}

pub struct BondHarmonic {
    k: f32,
    r0: f32,
    atoms: [usize; 2],
}

impl BondHarmonic {
    pub fn new(k: f32, r0: f32, atoms: [usize; 2]) -> Box<BondHarmonic> {
        Box::new(BondHarmonic { k, r0, atoms })
    }
}

impl BondedInteraction for BondHarmonic {
    fn calc(&self, ri: &Rvec, rj: &Rvec) -> (f32, [Rvec; 2]) {
        functions::bond_harm(self.k, self.r0, ri, rj)
    }
    fn atoms(&self) -> [usize; 2] {
        self.atoms
    }
}

pub struct AngleHarmonic {
    k: f32,
    t0: f32,
    atoms: [usize; 3],
}

impl AngleHarmonic {
    pub fn new(k: f32, t0: f32, atoms: [usize; 3]) -> Box<AngleHarmonic> {
        Box::new(AngleHarmonic { k, t0, atoms })
    }
}

impl AngleInteraction for AngleHarmonic {
    fn calc(&self, ri: &Rvec, rj: &Rvec, rk: &Rvec) -> (f32, [Rvec; 3]) {
        functions::angle_harm(self.k, self.t0, ri, rj, rk)
    }
    fn atoms(&self) -> [usize; 3] {
        self.atoms
    }
}

pub struct DihedralPeriodic {
    k: f32,
    n: f32,
    p0: f32,
    atoms: [usize; 4],
}

impl DihedralPeriodic {
    pub fn new(k: f32, n: f32, p0: f32, atoms: [usize; 4]) -> Box<DihedralPeriodic> {
        Box::new(DihedralPeriodic { k, n, p0, atoms })
    }
}

impl TorsionInteraction for DihedralPeriodic {
    fn calc(&self, ri: &Rvec, rj: &Rvec, rk: &Rvec, rl: &Rvec) -> (f32, [Rvec; 4]) {
        functions::pdih(self.k, self.n, self.p0, ri, rj, rk, rl)
    }
    fn atoms(&self) -> [usize; 4] {
        self.atoms
    }
}

pub struct ImproperDihedralHarmonic {
    k: f32,
    p0: f32,
    atoms: [usize; 4],
}

impl ImproperDihedralHarmonic {
    pub fn new(k: f32, p0: f32, atoms: [usize; 4]) -> Box<ImproperDihedralHarmonic> {
        Box::new(ImproperDihedralHarmonic { k, p0, atoms })
    }
}

impl TorsionInteraction for ImproperDihedralHarmonic {
    fn calc(&self, ri: &Rvec, rj: &Rvec, rk: &Rvec, rl: &Rvec) -> (f32, [Rvec; 4]) {
        functions::idih_harm(self.k, self.p0, ri, rj, rk, rl)
    }
    fn atoms(&self) -> [usize; 4] {
        self.atoms
    }
}

#[cfg(test)]
mod tests {
    use super::*;
}
