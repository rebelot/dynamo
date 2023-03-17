use crate::{topology::Topology, ffield::Forces};
mod reader;

#[derive(Debug)]
pub struct System {
    pub topology: Topology,
    pub forces: Forces,
    pub pbc: [f64; 3],
    pub potential: f64,
    pub press: f64,
}

impl System {
    pub fn new(topology: Topology, forces: Forces, pbc: [f64; 3]) -> System {
        System {
            topology,
            forces,
            pbc,
            potential: 0.0,
            press: 0.0,
        }
    }
    pub fn read(filename: &str) -> Self {
        return reader::read(filename);
    }
}
