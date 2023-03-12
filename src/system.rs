use crate::topology::Topology;

pub struct System {
    pub topology: Topology,
    pub r#box: [f64; 3],
    pub potential: f64,
    pub temp: f64,
    pub press: f64,
}
