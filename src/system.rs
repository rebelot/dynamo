use crate::{Rvec, ffield::Forces, topology::Topology};
mod reader;

pub struct System {
    pub topology: Topology,
    pub forces: Forces,
    pub pbc: Rvec,
    pub potential: f64,
    pub press: f64,
}

impl System {
    pub fn new(topology: Topology, forces: Forces, pbc: Rvec) -> System {
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
