use crate::{ffield::Forces, topology::Topology, Rvec};

pub struct System {
    pub topology: Topology,
    pub forces: Forces,
    pub pbc: Rvec,
    pub potential: f32,
    pub press: f32,
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
}
