use crate::ffield::FFList;
use std::{cell::RefCell, rc::Rc};

pub mod atom;
mod reader;

pub fn test() {
    println!("Hello from topology!");
}

#[derive(Debug)]
pub struct Topology {
    pub atoms: Vec<Rc<RefCell<atom::Atom>>>,
    pub ff: FFList,
}

impl Topology {
    pub fn new() -> Topology {
        Topology {
            atoms: Vec::new(),
            ff: FFList::new(),
        }
    }
    pub fn add_bond(&mut self, a1: usize, a2: usize) {
        self.atoms[a1]
            .borrow_mut()
            .bonds
            .push(Rc::clone(&self.atoms[a2]));
        self.atoms[a2]
            .borrow_mut()
            .bonds
            .push(Rc::clone(&self.atoms[a1]));
    }
    pub fn add_atom(
        &mut self,
        atomtype: String,
        name: String,
        element: String,
        mass: f64,
        vdw: f64,
        charge: f64,
        coords: [f64; 3],
    ) -> usize {
        let index = self.atoms.len();
        let atom = atom::Atom::new(index, atomtype, name, element, mass, vdw, charge, coords);
        self.atoms.push(atom);
        return index;
    }
    pub fn read(filename: &str) -> Topology {
        return reader::read(filename);
    }
}
