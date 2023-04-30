use rayon::prelude::*;
use std::sync::{mpsc, Arc, Mutex};
use std::thread;

#[derive(Debug, Clone)]
struct Atom {
    pos: [f32; 3],
    f: [f32; 3],
}

impl Atom {
    fn new() -> Atom {
        Atom {
            pos: [0.0, 0.0, 0.0],
            f: [0.0, 0.0, 0.0],
        }
    }
}

trait Bond {
    fn atoms(&self) -> [usize; 2];
    fn calc(&self, ri: &[f32; 3], rj: &[f32; 3]) -> (f32, [[f32; 3]; 2]);
}

struct BType1 {
    atoms: [usize; 2],
}
impl BType1 {
    fn new(i: usize, j: usize) -> BType1 {
        BType1 { atoms: [i, j] }
    }
}

impl Bond for BType1 {
    fn calc(&self, ri: &[f32; 3], rj: &[f32; 3]) -> (f32, [[f32; 3]; 2]) {
        (1.0, [[1.0, 1.0, 1.0], [-1.0, -1.0, -1.0]])
    }
    fn atoms(&self) -> [usize; 2] {
        self.atoms
    }
}

struct BType2 {
    atoms: [usize; 2],
}
impl BType2 {
    fn new(i: usize, j: usize) -> BType2 {
        BType2 { atoms: [i, j] }
    }
}

impl Bond for BType2 {
    fn calc(&self, ri: &[f32; 3], rj: &[f32; 3]) -> (f32, [[f32; 3]; 2]) {
        (1.0, [[1.0, 1.0, 1.0], [-1.0, -1.0, -1.0]])
    }
    fn atoms(&self) -> [usize; 2] {
        self.atoms
    }
}

fn main() {
    let mut atoms = Vec::new();
    atoms.push(Mutex::new(Atom::new()));
    atoms.push(Mutex::new(Atom::new()));
    atoms.push(Mutex::new(Atom::new()));

    let mut bonds = Vec::<Box<dyn Bond + Send + Sync>>::new();
    bonds.push(Box::new(BType1::new(0, 1)));
    bonds.push(Box::new(BType2::new(1, 2)));

    bonds.par_iter().for_each(|bond| {
        let [i, j] = bond.atoms();
        let mut ai = atoms[i].lock().unwrap();
        let mut aj = atoms[j].lock().unwrap();

        let (u, f) = bond.calc(&ai.pos, &aj.pos);

        for i in 0..3 {
            ai.f[i] += f[0][i];
            aj.f[i] += f[1][i];
        }
    });

    println!("{:?}", atoms);
}
