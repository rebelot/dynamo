use crate::ffield;

#[derive(Debug)]
pub struct Atom {
    pub index: usize,
    pub atomtype: String,
    pub name: String,
    pub element: String,
    pub pos: [f64; 3],
    pub vel: [f64; 3],
    pub force: [f64; 3],
    pub prev_force: [f64; 3],
    pub mass: f64,
    pub vdw: f64,
    pub charge: f64,
    pub lj: ffield::LJParams,
}

impl Atom {
    pub fn new(
        index: usize,
        atomtype: String,
        name: String,
        element: String,
        mass: f64,
        vdw: f64,
        charge: f64,
        pos: [f64; 3],
    ) -> Atom {
        Atom {
            index,
            atomtype,
            name,
            element,
            mass,
            vdw,
            charge,
            pos,
            vel: [0.0, 0.0, 0.0],
            force: [0.0, 0.0, 0.0],
            prev_force: [0.0, 0.0, 0.0],
            lj: ffield::LJParams::new(0.0, 0.0),
        }
    }
}

