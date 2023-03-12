use std::rc::Rc;

use crate::{system::System, topology::atom::Atom};

const DIM: usize = 3;

pub struct VelocityVerlet {
    pub system: System,
    pub dt: f64,
    pub time: f64,
}

impl VelocityVerlet {
    pub fn new(system: System, dt: f64) -> VelocityVerlet {
        VelocityVerlet { system, dt, time: 0.0 }
    }
    pub fn step(&mut self) {
        self.time += self.dt;

        for a in self.system.topology.atoms.iter() {
            update_pos(&mut a.borrow_mut(), self.dt);
        }

        self.system.topology.ff.calc();

        for a in self.system.topology.atoms.iter() {
            update_vel(&mut a.borrow_mut(), self.dt);
        }

    }
}

fn update_pos(a: &mut Atom, dt: f64){
    for i in 0..DIM {
        a.pos[i] += a.vel[i] * dt + 0.5 * a.force[i] / a.mass * dt * dt;
    }
}

fn update_vel(a: &mut Atom, dt: f64){
    for i in 0..DIM {
        a.vel[i] += 0.5 * (a.force[i] + a.force[i]) / a.mass * dt;
    }

}

