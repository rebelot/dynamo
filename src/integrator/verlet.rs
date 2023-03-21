use crate::{DIM, system::System, topology::atom::Atom};

pub struct VelocityVerlet {
    pub dt: f32,
    pub time: f32,
    pub step: i32,
}

impl VelocityVerlet {
    pub fn new(dt: f32) -> VelocityVerlet {
        VelocityVerlet {
            dt,
            time: 0.0,
            step: 0,
        }
    }
    pub fn step(&mut self, sys: &mut System) {
        self.time += self.dt;

        for a in sys.topology.atoms.iter_mut() {
            Self::update_pos(&self.dt, a);
        }

        sys.forces.calc(&mut sys.topology.atoms);

        for a in sys.topology.atoms.iter_mut() {
            Self::update_vel(&self.dt, a);
        }

        sys.topology.init_forces();
    }

    fn update_pos(dt: &f32, a: &mut Atom) {
        for i in 0..DIM {
            a.pos[i] += a.vel[i] * dt + 0.5 * a.prev_force[i] / a.mass * dt * dt;
        }
    }

    fn update_vel(dt: &f32, a: &mut Atom) {
        for i in 0..DIM {
            a.vel[i] += 0.5 * (a.force[i] + a.prev_force[i]) / a.mass * dt;
        }
    }
}
