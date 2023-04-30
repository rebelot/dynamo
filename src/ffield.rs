// use crate::topology::atom::Atom;
use crate::{Rvec, DIM};
// use rayon::prelude::*;

mod functions;

pub struct Forces {
    pub bonds: Vec<Box<dyn BondedInteraction + Send + Sync>>,
    pub angles: Vec<Box<dyn AngleInteraction + Send + Sync>>,
    pub torsions: Vec<Box<dyn TorsionInteraction + Send + Sync>>,
    pub pairs: Vec<Box<dyn BondedInteraction + Send + Sync>>,
}

impl Default for Forces {
    fn default() -> Self {
        Forces::new()
    }
}

impl Forces {
    pub fn new() -> Forces {
        Forces {
            bonds: Vec::new(),
            angles: Vec::new(),
            torsions: Vec::new(),
            pairs: Vec::new(),
        }
    }
    pub fn parse_interaction(&mut self, interaction: (String, Vec<f32>), offset: usize) {
        match interaction.0.as_str() {
            "bond_harm" => {
                let mut atoms: [usize; 2] = [0, 0];
                interaction.1[..2]
                    .iter()
                    .enumerate()
                    .for_each(|(i, x)| atoms[i] = *x as usize + offset);
                let k = interaction.1[2];
                let r0 = interaction.1[3];
                self.bonds.push(BondHarmonic::new(k, r0, atoms));
            }
            "angle_harm" => {
                let mut atoms: [usize; 3] = [0, 0, 0];
                interaction.1[..3]
                    .iter()
                    .enumerate()
                    .for_each(|(i, x)| atoms[i] = *x as usize + offset);
                let k = interaction.1[3];
                let t0 = interaction.1[4];
                self.angles.push(AngleHarmonic::new(k, t0, atoms));
            }

            "pdih" => {
                let mut atoms: [usize; 4] = [0, 0, 0, 0];
                interaction.1[..4]
                    .iter()
                    .enumerate()
                    .for_each(|(i, x)| atoms[i] = *x as usize + offset);
                let k = interaction.1[4];
                let n = interaction.1[5];
                let p0 = interaction.1[6];
                self.torsions.push(DihedralPeriodic::new(k, n, p0, atoms));
            }

            "idih_harm" => {
                let mut atoms: [usize; 4] = [0, 0, 0, 0];
                interaction.1[..4]
                    .iter()
                    .enumerate()
                    .for_each(|(i, x)| atoms[i] = *x as usize + offset);
                let k = interaction.1[4];
                let p0 = interaction.1[5];
                self.torsions
                    .push(ImproperDihedralHarmonic::new(k, p0, atoms));
            }
            "lj_pair" => {
                let mut atoms: [usize; 2] = [0, 0];
                interaction.1[..2]
                    .iter()
                    .enumerate()
                    .for_each(|(i, x)| atoms[i] = *x as usize + offset);
                let v = interaction.1[2];
                let w = interaction.1[3];
                self.pairs.push(BondedLJPair::new(v, w, atoms));
            }
            _ => panic!("Unknown interaction"),
        }
    }

    pub fn calc(&self, positions: &[Rvec], forces: &mut [Rvec]) -> f32 {
        let mut tot_u = 0.0;

        self.bonds.iter().for_each(|bond| {
            let [i, j] = bond.atoms();
            let (u, f) = bond.calc(&positions[i], &positions[j]);

            // tot_u += u;
            for d in 0..DIM {
                forces[i][d] += f[0][d];
                forces[j][d] += f[1][d];
            }
        });

        self.angles.iter().for_each(|angle| {
            let [i, j, k] = angle.atoms();
            let (u, f) = angle.calc(&positions[i], &positions[j], &positions[k]);

            tot_u += u;
            for d in 0..DIM {
                forces[i][d] += f[0][d];
                forces[j][d] += f[1][d];
                forces[k][d] += f[2][d];
            }
        });
        self.torsions.iter().for_each(|torsion| {
            let [i, j, k, l] = torsion.atoms();
            let (u, f) = torsion.calc(&positions[i], &positions[j], &positions[k], &positions[l]);

            tot_u += u;
            for d in 0..DIM {
                forces[i][d] += f[0][d];
                forces[j][d] += f[1][d];
                forces[k][d] += f[2][d];
                forces[l][d] += f[3][d];
            }
        });
        tot_u
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

pub struct BondedLJPair {
    v: f32,
    w: f32,
    atoms: [usize; 2],
}

impl BondedLJPair {
    pub fn new(v: f32, w: f32, atoms: [usize; 2]) -> Box<BondedLJPair> {
        Box::new(BondedLJPair { v, w, atoms })
    }
}

impl BondedInteraction for BondedLJPair {
    fn calc(&self, ri: &Rvec, rj: &Rvec) -> (f32, [Rvec; 2]) {
        functions::lj(self.v, self.w, ri, rj)
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
