use crate::ffield;
use crate::Rvec;

#[derive(Debug)]
pub struct Atom {
    pub index: usize,
    pub atomtype: String,
    pub name: String,
    pub element: String,
    pub pos: Rvec,
    pub vel: Rvec,
    pub force: Rvec,
    pub prev_force: Rvec,
    pub mass: f32,
    pub vdw: f32,
    pub charge: f32,
    // pub lj: ffield::LJParams,
}

impl Atom {
    pub fn new(
        index: usize,
        atomtype: String,
        name: String,
        element: String,
        mass: f32,
        vdw: f32,
        charge: f32,
        pos: Rvec,
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
            // lj: ffield::LJParams::new(0.0, 0.0),
        }
    }
}
