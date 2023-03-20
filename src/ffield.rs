use crate::topology::atom::Atom;

mod functions;

const DIM: usize = 3;

pub struct Forces {
    pub bonds: Vec<Box<dyn Interaction>>,
    pub angles: Vec<Box<dyn Interaction>>,
    pub torsions: Vec<Box<dyn Interaction>>,
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
        for b in self.bonds.iter() {
            b.calc(atoms);
        }
        for a in self.angles.iter() {
            a.calc(atoms);
        }
        for t in self.torsions.iter() {
            t.calc(atoms);
        }
    }
}

pub trait Interaction {
    fn calc(&self, atoms: &mut Vec<Atom>);
}

pub struct BondHarmonic {
    k: f64,
    r0: f64,
    atoms: [usize; 2],
}

impl BondHarmonic {
    pub fn new(k: f64, r0: f64, atoms: [usize; 2]) -> Box<BondHarmonic> {
        Box::new(BondHarmonic { k, r0, atoms })
    }
}

impl Interaction for BondHarmonic {
    fn calc(&self, atoms: &mut Vec<Atom>) {
        let r = [&atoms[self.atoms[0]].pos, &atoms[self.atoms[1]].pos];

        let (u, f) = functions::bond_harm(&self.k, &self.r0, r);

        for i in 0..DIM {
            atoms[self.atoms[0]].force[i] += f[0][i];
            atoms[self.atoms[1]].force[i] += f[1][i];
        }
    }
}

pub struct AngleHarmonic {
    k: f64,
    t0: f64,
    atoms: [usize; 3],
}

impl AngleHarmonic {
    pub fn new(k: f64, t0: f64, atoms: [usize; 3]) -> Box<AngleHarmonic> {
        Box::new(AngleHarmonic { k, t0, atoms })
    }
}

impl Interaction for AngleHarmonic {
    fn calc(&self, atoms: &mut Vec<Atom>) {
        let r = [
            &atoms[self.atoms[0]].pos,
            &atoms[self.atoms[1]].pos,
            &atoms[self.atoms[2]].pos,
        ];

        let (u, f) = functions::angle_harm(&self.k, &self.t0, r);

        for i in 0..DIM {
            atoms[self.atoms[0]].force[i] += f[0][i];
            atoms[self.atoms[1]].force[i] += f[1][i];
            atoms[self.atoms[2]].force[i] += f[2][i];
        }
    }
}

pub struct DihedralPeriodic {
    k: f64,
    n: f64,
    p0: f64,
    atoms: [usize; 4],
}

impl DihedralPeriodic {
    pub fn new(k: f64, n: f64, p0: f64, atoms: [usize; 4]) -> Box<DihedralPeriodic> {
        Box::new(DihedralPeriodic { k, n, p0, atoms })
    }
}

impl Interaction for DihedralPeriodic {
    fn calc(&self, atoms: &mut Vec<Atom>) {
        let r = [
            &atoms[self.atoms[0]].pos,
            &atoms[self.atoms[1]].pos,
            &atoms[self.atoms[2]].pos,
            &atoms[self.atoms[3]].pos,
        ];

        let (u, f) = functions::pdih(&self.k, &self.n, &self.p0, r);

        for i in 0..DIM {
            atoms[self.atoms[0]].force[i] += f[0][i];
            atoms[self.atoms[1]].force[i] += f[1][i];
            atoms[self.atoms[2]].force[i] += f[2][i];
            atoms[self.atoms[2]].force[i] += f[2][i];
        }
    }
}

pub struct ImproperDihedralHarmonic {
    k: f64,
    p0: f64,
    atoms: [usize; 4],
}

impl ImproperDihedralHarmonic {
    pub fn new(k: f64, p0: f64, atoms: [usize; 4]) -> Box<ImproperDihedralHarmonic> {
        Box::new(ImproperDihedralHarmonic { k, p0, atoms })
    }
}

impl Interaction for ImproperDihedralHarmonic {
    fn calc(&self, atoms: &mut Vec<Atom>) {
        let r = [
            &atoms[self.atoms[0]].pos,
            &atoms[self.atoms[1]].pos,
            &atoms[self.atoms[2]].pos,
            &atoms[self.atoms[3]].pos,
        ];

        let (u, f) = functions::idih_harm(&self.k, &self.p0, r);

        for i in 0..DIM {
            atoms[self.atoms[0]].force[i] += f[0][i];
            atoms[self.atoms[1]].force[i] += f[1][i];
            atoms[self.atoms[2]].force[i] += f[2][i];
            atoms[self.atoms[2]].force[i] += f[2][i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
}
