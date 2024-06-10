use crate::{ffield::Forces, topology::Topology, Rvec};

pub struct System {
    pub topology: Topology,
    pub forces: Forces,
    pub potential: f32,
    pub press: f32,
    pub temp: f32,
}

impl System {
    pub fn new(topology: Topology, forces: Forces) -> System {
        System {
            topology,
            forces,
            potential: 0.0,
            press: 0.0,
            temp: 0.0,
        }
    }
}
