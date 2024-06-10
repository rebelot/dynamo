use crate::linalg::DEG2RAD;
use std::ops::Index;

// use crate::topology::atom::Atom;
use crate::{
    topology::{atom::Atom, molecule::Molecule, Topology},
    Rvec, DIM,
};
// use rayon::prelude::*;

mod functions;

pub struct Forces {
    pub bonds: Vec<Box<dyn TwoAtomInteraction + Send + Sync>>,
    pub angles: Vec<Box<dyn ThreeAtomInteraction + Send + Sync>>,
    pub torsions: Vec<Box<dyn FourAtomInteraction + Send + Sync>>,
    pub pairs: Vec<Box<dyn TwoAtomInteraction + Send + Sync>>,
    comb_rule: fn(f32, f32, f32, f32) -> (f32, f32),
    qqscale: f32,
    ljscale: f32,
}

impl Forces {
    #[allow(clippy::unnecessary_operation)]
    pub fn new(top: &Topology) -> Forces {
        let mut ff = Forces {
            bonds: Vec::new(),
            angles: Vec::new(),
            torsions: Vec::new(),
            pairs: Vec::new(),
            comb_rule: match top.defaults.comb_rule.as_str() {
                "geom" => functions::comb_rule_geom,
                "LB" => functions::comb_rule_LB,
                wtf => panic!("Unknown combination rule: {}", wtf),
            },
            qqscale: top.defaults.qqscale.unwrap_or(1.0),
            ljscale: top.defaults.ljscale.unwrap_or(1.0),
        };
        ff.build(top);
        ff
    }

    fn build(&mut self, top: &Topology) {
        let mut off = 0;
        for mol in &top.molecules {
            let natoms = mol.atoms.len();
            for _ in 0..mol.nmols {
                for interaction in &mol.bonded_interactions {
                    self.parse_interaction(mol, interaction, off);
                }
                off += natoms;
            }
        }
    }

    fn parse_interaction(&mut self, mol: &Molecule, interaction: &str, offset: usize) {
        let fields = interaction.split_whitespace().collect::<Vec<&str>>();
        let funct = fields[0];
        let params = &fields[1..];
        match funct {
            "bond_harm" => {
                let mut atoms: [usize; 2] = [0, 0];
                params[..2]
                    .iter()
                    .enumerate()
                    .for_each(|(i, x)| atoms[i] = x.parse::<usize>().unwrap() + offset - 1);
                let r0 = params[2].parse::<f32>().unwrap();
                let k = params[3].parse::<f32>().unwrap();
                self.bonds.push(BondHarmonic::new(k, r0, atoms));
            }

            "angle_harm" => {
                let mut atoms: [usize; 3] = [0, 0, 0];
                params[..3]
                    .iter()
                    .enumerate()
                    .for_each(|(i, x)| atoms[i] = x.parse::<usize>().unwrap() + offset - 1);
                let t0 = params[3].parse::<f32>().unwrap() * DEG2RAD;
                let k = params[4].parse::<f32>().unwrap();
                self.angles.push(AngleHarmonic::new(k, t0, atoms));
            }

            "pdih" => {
                let mut atoms: [usize; 4] = [0, 0, 0, 0];
                params[..4]
                    .iter()
                    .enumerate()
                    .for_each(|(i, x)| atoms[i] = x.parse::<usize>().unwrap() + offset - 1);
                let p0 = params[4].parse::<f32>().unwrap() * DEG2RAD;
                let k = params[5].parse::<f32>().unwrap();
                let n = params[6].parse::<f32>().unwrap();
                self.torsions.push(DihedralPeriodic::new(k, n, p0, atoms));
            }

            "idih_harm" => {
                let mut atoms: [usize; 4] = [0, 0, 0, 0];
                params[..4]
                    .iter()
                    .enumerate()
                    .for_each(|(i, x)| atoms[i] = x.parse::<usize>().unwrap() + offset - 1);
                let p0 = params[4].parse::<f32>().unwrap() * DEG2RAD;
                let k = params[5].parse::<f32>().unwrap();
                self.torsions
                    .push(ImproperDihedralHarmonic::new(k, p0, atoms));
            }
            "lj_pair" => {
                let mut atoms: [usize; 2] = [0, 0];
                params[..2]
                    .iter()
                    .enumerate()
                    .for_each(|(i, x)| atoms[i] = x.parse::<usize>().unwrap() + offset - 1);
                let v = params.get(2);
                let w = params.get(3);

                if let (Some(v), Some(w)) = (v, w) {
                    self.pairs.push(BondedLJPair::new(
                        v.parse().unwrap(),
                        w.parse().unwrap(),
                        atoms,
                    ));
                } else {
                    let [mut ai, mut aj] = atoms;
                    ai -= offset;
                    aj -= offset;
                    let (vcomb, wcomb) = (self.comb_rule)(
                        mol.atoms[ai].v,
                        mol.atoms[ai].w,
                        mol.atoms[aj].v,
                        mol.atoms[aj].w,
                    );
                    self.pairs.push(BondedLJPair::new(vcomb, wcomb, atoms));
                }
            }
            wtf => panic!("Unknown interaction: {}", wtf),
        }
    }

    pub fn calc(&self, positions: &[Rvec], forces: &mut [Rvec]) -> f32 {
        let mut tot_u = 0.0;

        self.bonds.iter().for_each(|bond| {
            let [i, j] = bond.atoms();
            let (u, f) = bond.calc(&positions[i], &positions[j]);

            tot_u += u;
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
        self.pairs.iter().for_each(|pair| {
            let [i, j] = pair.atoms();
            let (u, f) = pair.calc(&positions[i], &positions[j]);

            tot_u += u * self.ljscale;
            for d in 0..DIM {
                forces[i][d] += f[0][d] * self.ljscale;
                forces[j][d] += f[1][d] * self.ljscale;
            }
        });
        tot_u
    }
}

pub trait TwoAtomInteraction {
    fn calc(&self, ri: &Rvec, rj: &Rvec) -> (f32, [Rvec; 2]);
    fn atoms(&self) -> [usize; 2];
}

pub trait ThreeAtomInteraction {
    fn calc(&self, ri: &Rvec, rj: &Rvec, rk: &Rvec) -> (f32, [Rvec; 3]);
    fn atoms(&self) -> [usize; 3];
}

pub trait FourAtomInteraction {
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

impl TwoAtomInteraction for BondHarmonic {
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

impl TwoAtomInteraction for BondedLJPair {
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

impl ThreeAtomInteraction for AngleHarmonic {
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

impl FourAtomInteraction for DihedralPeriodic {
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

impl FourAtomInteraction for ImproperDihedralHarmonic {
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

    #[test]
    fn it_builds() {
        //let top = Topology::read("tests/diala.top");
        //let ff = Forces::new(&top);
        let mut top = Topology::new();
        top.defaults.comb_rule = "geom".to_string();
        top.add_atomtype("a".to_string(), 0, 0.0, 0.0, 0.0);

        top.add_molecule("one".to_string(), 1, 3);
        top.add_atom(0, "a", "a", 1, "a", 0.0);
        top.add_atom(0, "a", "a", 1, "a", 0.0);
        top.add_atom(0, "a", "a", 1, "a", 0.0);
        top.add_bonded_interaction(0, "bond_harm 1 2 0 0");
        top.add_bonded_interaction(0, "bond_harm 1 3 0 0");
        top.add_bonded_interaction(0, "bond_harm 2 3 0 0");

        top.add_molecule("two".to_string(), 10, 3);
        top.add_atom(1, "a", "a", 1, "a", 0.0);
        top.add_atom(1, "a", "a", 1, "a", 0.0);
        top.add_atom(1, "a", "a", 1, "a", 0.0);
        top.add_bonded_interaction(1, "bond_harm 1 2 0 0");
        top.add_bonded_interaction(1, "bond_harm 1 3 0 0");
        top.add_bonded_interaction(1, "bond_harm 2 3 0 0");

        let ff = Forces::new(&top);
        let [a1, a2] = ff.bonds[ff.bonds.len() - 1].atoms();
        // 3 + 3*10 atoms = 33; (2, 3) is (32, 33), indexed at (31, 32)
        assert_eq!(a1, 31);
        assert_eq!(a2, 32);
    }
}
