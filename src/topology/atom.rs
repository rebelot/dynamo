use std::cell::RefCell;
use std::fmt::Debug;
use std::rc::Rc;

use crate::ffield;

pub struct Atom {
    pub index: usize,
    pub atomtype: String,
    pub name: String,
    pub element: String,
    pub pos: [f64; 3],
    pub vel: [f64; 3],
    pub force: [f64; 3],
    pub mass: f64,
    pub vdw: f64,
    pub charge: f64,
    pub LJ: ffield::LJParams,
    pub bonds: Vec<Rc<RefCell<Atom>>>,
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
    ) -> Rc<RefCell<Atom>> {
        Rc::new(RefCell::new(Atom {
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
            LJ: ffield::LJParams::new(0.0, 0.0),
            bonds: Vec::new(),
        }))
    }
}

impl Debug for Atom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let bonds: Vec<usize> = self.bonds.iter().map(|x| x.borrow().index).collect();
        f.debug_struct("Atom")
            .field("index", &self.index)
            .field("atomtype", &self.atomtype)
            .field("name", &self.name)
            .field("element", &self.element)
            .field("pos", &self.pos)
            .field("mass", &self.mass)
            .field("vdw", &self.vdw)
            .field("charge", &self.charge)
            .field("vel", &self.vel)
            .field("force", &self.force)
            .field("LJ", &self.LJ)
            .field("bonds", &bonds)
            .finish()
    }
}
