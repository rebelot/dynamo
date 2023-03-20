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
        Self { k, r0, atoms }
    }

    fn calc(&mut self, atoms: &mut Vec<Atom>) {
        let ri = &atoms[self.atoms[0]].pos;
        let rj = &atoms[self.atoms[1]].pos;

        let rij = displace_vec(ri, rj);
        let norm2_rij = norm2(&rij);
        let norm_rij = norm2_rij.sqrt();

        let (u, mut f) = harmonic(&self.k, &self.r0, &norm_rij);
        f *= norm_rij.recip();
        // *f *= norm2_rij.sqrt().recip(); // will this be faster?

        let mut dri;
        for i in 0..DIM {
            dri = f * -rij[i];
            atoms[self.atoms[0]].force[i] += dri;
            atoms[self.atoms[1]].force[i] -= dri;
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
        Self { k, t0, atoms }
    }
    fn calc(&mut self, atoms: &mut Vec<Atom>) {
        let ri = &atoms[self.atoms[0]].pos;
        let rj = &atoms[self.atoms[1]].pos;
        let rk = &atoms[self.atoms[2]].pos;

        let rji = displace_vec(rj, ri);
        let rjk = displace_vec(rj, rk);

        let norm_rji_rjk_1 = (norm2(&rji) * norm2(&rjk)).sqrt().recip();

        let cos_t = dot(&rji, &rjk) * norm_rji_rjk_1;
        let t = cos_t.acos();

        let (u, mut f) = harmonic(&self.k, &self.t0, &t);

        f *= -(1.0 - cos_t * cos_t).sqrt().recip() * norm_rji_rjk_1;
        let (mut dri, mut drj, mut drk);
        for i in 0..DIM {
            dri = f * rjk[i];
            drk = f * rji[i];
            drj = -dri - drk;
            atoms[self.atoms[0]].force[i] += dri;
            atoms[self.atoms[1]].force[i] += drj;
            atoms[self.atoms[2]].force[i] += drk;
        }
    }
}

#[derive(Debug)]
pub struct Dihedral {
    k: f64,
    n: f64,
    t0: f64,
    atoms: [usize; 4],
}

impl Dihedral {
    pub fn new(k: f64, t0: f64, n: f64, atoms: [usize; 4]) -> Self {
        Self { k, t0, n, atoms }
    }
    fn calc(&mut self, atoms: &mut Vec<Atom>) {
        let psi = self.partial(atoms);
        let (u, f) = if self.k != 0.0 {
            periodic(&self.k, &self.t0, &psi, &self.n)
        } else {
            harmonic(&self.k, &self.t0, &psi)
        };
        for i in 0..DIM {
            atoms[self.atoms[0]].force[i] += f;
            atoms[self.atoms[1]].force[i] += f;
            atoms[self.atoms[2]].force[i] += f;
            atoms[self.atoms[3]].force[i] += f;
        }
    }
    fn partial(&mut self, atoms: &mut Vec<Atom>) -> f64 {
        let ri = &atoms[self.atoms[0]].pos;
        let rj = &atoms[self.atoms[1]].pos;
        let rk = &atoms[self.atoms[2]].pos;
        let rl = &atoms[self.atoms[3]].pos;

        let rij = displace_vec(ri, rj);
        let rjk = displace_vec(rj, rk);
        let rkl = displace_vec(rk, rl);

        let nijk = cross(&rij, &rjk);
        let njkl = cross(&rjk, &rkl);

        let y = dot(&cross(&nijk, &njkl), &rjk) * norm2(&rjk).sqrt().recip();
        let x = dot(&nijk, &njkl);
        let psi = y.atan2(x);

        let f = (1.0 + (y / x).powi(2)) * y * x.sqrt().recip();

        let dnijk_drij = [rjk[2] - rjk[1], rjk[0] - rjk[2], rjk[1] - rjk[0]];
        let dnijk_drjk = [rij[1] - rij[2], rij[2] - rij[0], rij[0] - rij[1]];
        let dnjkl_drjk = [rkl[2] - rkl[1], rkl[0] - rkl[2], rkl[1] - rkl[0]];
        let dnjkl_drkl = [rjk[1] - rjk[2], rjk[2] - rjk[0], rjk[0] - rjk[1]];

        let (mut drij, mut drjk, mut drkl);
        for i in 0..DIM {
            drij = f * njkl[i] * dnijk_drij[i];
            drjk = f * (njkl[i] * dnijk_drjk[i] + nijk[i] * dnjkl_drjk[i]);
            drkl = f * nijk[i] * dnjkl_drkl[i];
            atoms[self.atoms[0]].force[i] -= drij;
            atoms[self.atoms[1]].force[i] += drij - drjk;
            atoms[self.atoms[2]].force[i] += drjk - drkl;
            atoms[self.atoms[3]].force[i] += drkl;
        }
        return psi;
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

fn harmonic(k: &f64, x0: &f64, x: &f64) -> (f64, f64) {
    let dx = x - x0;
    let kdx = k * dx;
    return (0.5 * kdx * dx, -kdx);
}

fn periodic(k: &f64, x0: &f64, x: &f64, n: &f64) -> (f64, f64) {
    let dx = x * n - x0;
    return (k * (1.0 + dx.cos()), k * n * dx.sin());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_harmonic() {
        let k = 1.0;
        let r0 = 0.0;
        let r = 1.0;
        let (u, f) = harmonic(&k, &r0, &r);
        assert_eq!(u, 0.5);
        assert_eq!(f, -1.0);
    }

    #[test]
    fn test_periodic() {
        let k = 1.0;
        let r0 = 0.0;
        let r = 90_f64.to_radians();
        let n = 1.0;
        let (u, f) = periodic(&k, &r0, &r, &n);
        assert_eq!(u, 1.0);
        assert_eq!(f, 1.0);
    }

    #[test]
    fn test_bond() {
        let mut atoms = Vec::<Atom>::new();
        atoms.push(Atom::new(
            0,
            String::from(""),
            String::from(""),
            String::from(""),
            1.0,
            0.0,
            0.0,
            [0.0, 0.0, 0.0],
        ));
        atoms.push(Atom::new(
            1,
            String::from(""),
            String::from(""),
            String::from(""),
            1.0,
            0.0,
            0.0,
            [2.0, 0.0, 0.0],
        ));

        let mut forces = BondedForces::new();
        forces.bonds.push(Bond::new(1.0, 1.0, [0, 1]));

        forces.calc(&mut atoms);

        assert_eq!(
            atoms[0]
                .force
                .iter()
                .map(|x| format!("{:.3}", x))
                .collect::<Vec<String>>(),
            [
                String::from("1.000"),
                String::from("0.000"),
                String::from("0.000")
            ]
        );
        assert_eq!(
            atoms[1]
                .force
                .iter()
                .map(|x| format!("{:.3}", x))
                .collect::<Vec<String>>(),
            [
                String::from("-1.000"),
                String::from("0.000"),
                String::from("0.000")
            ]
        );
    }

    #[test]
    fn test_angle() {
        let mut atoms = Vec::<Atom>::new();
        atoms.push(Atom::new(
            0,
            String::from(""),
            String::from(""),
            String::from(""),
            1.0,
            0.0,
            0.0,
            [0.0, 1.0, 0.0],
        ));
        atoms.push(Atom::new(
            1,
            String::from(""),
            String::from(""),
            String::from(""),
            1.0,
            0.0,
            0.0,
            [0.0, 0.0, 0.0],
        ));

        atoms.push(Atom::new(
            2,
            String::from(""),
            String::from(""),
            String::from(""),
            1.0,
            0.0,
            0.0,
            [1.0, 0.0, 0.0],
        ));

        let mut forces = BondedForces::new();
        forces
            .angles
            .push(Angle::new(1.0, 0_f64.to_radians(), [0, 1, 2]));

        forces.calc(&mut atoms);

        println!("{:?}", atoms[0].force);
        println!("{:?}", atoms[1].force);
        println!("{:?}", atoms[2].force);
        assert_eq!(
            atoms[0]
                .force
                .iter()
                .map(|x| format!("{:.3}", x))
                .collect::<Vec<String>>(),
            [
                String::from("1.571"),
                String::from("0.000"),
                String::from("0.000")
            ]
        );
        assert_eq!(
            atoms[1]
                .force
                .iter()
                .map(|x| format!("{:.3}", x))
                .collect::<Vec<String>>(),
            [
                String::from("-1.571"),
                String::from("-1.571"),
                String::from("0.000")
            ]
        );
        assert_eq!(
            atoms[2]
                .force
                .iter()
                .map(|x| format!("{:.3}", x))
                .collect::<Vec<String>>(),
            [
                String::from("0.000"),
                String::from("1.571"),
                String::from("0.000")
            ]
        );
    }

    #[test]
    fn test_dihedral() {}
}
