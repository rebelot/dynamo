// use rayon::prelude::*;

use crate::{ffield::Forces, system::System, topology::atom::Atom, Rvec, DIM};

pub struct VelocityVerlet {
    pub dt: f32,
    pub time: f32,
    pub step: i32,
    cache: Vec<Rvec>,
}

impl VelocityVerlet {
    pub fn new(dt: f32, n: usize) -> VelocityVerlet {
        VelocityVerlet {
            dt,
            time: 0.0,
            step: 0,
            cache: vec![[0.0; DIM]; n],
        }
    }
    pub fn step(
        &mut self,
        ffield: &Forces,
        coords: &mut [Rvec],
        forces: &mut [Rvec],
        velocities: &mut [Rvec],
        masses: &[f32],
    ) -> f32 {
        self.time += self.dt;
        self.step += 1;

        coords
            .iter_mut()
            .zip(velocities.iter())
            .zip(forces.iter_mut())
            .zip(masses.iter())
            .for_each(|(((crd, vel), frc), mass)| Self::update_pos(&self.dt, crd, vel, frc, mass));

        self.cache.copy_from_slice(forces);
        forces.fill([0.0; DIM]);
        let u = ffield.calc(coords, forces);

        velocities
            .iter_mut()
            .zip(forces.iter())
            .zip(self.cache.iter())
            .zip(masses.iter())
            .for_each(|(((vel, frc), frcp), mass)| {
                Self::update_vel(&self.dt, vel, frc, frcp, mass)
            });
        u
    }

    #[inline]
    fn update_pos(dt: &f32, crd: &mut Rvec, vel: &Rvec, frc: &Rvec, mass: &f32) {
        for i in 0..DIM {
            crd[i] += vel[i] * dt + 0.5 * frc[i] / mass * dt * dt;
        }
    }

    #[inline]
    fn update_vel(dt: &f32, vel: &mut Rvec, frc: &Rvec, frcp: &Rvec, mass: &f32) {
        for i in 0..DIM {
            vel[i] += 0.5 * (frc[i] + frcp[i]) / mass * dt;
        }
    }
}
