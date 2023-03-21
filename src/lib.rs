pub mod ffield;
pub mod linalg;
pub mod system;
pub mod topology;
pub mod trajectory {
    pub mod writer;
}
pub mod integrator {
    pub mod verlet;
}

pub const DIM: usize = 3;
pub type Rvec = [f32; DIM];
