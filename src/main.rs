use dynamo::ffield::Forces;
use dynamo::integrator::verlet::VelocityVerlet;
use dynamo::system::System;
use dynamo::topology::Topology;
use dynamo::trajectory::writer::TrajectoryWriter;
use dynamo::{trajectory, Rvec};

fn main() {
    //rayon::ThreadPoolBuilder::new().num_threads(1).build_global().unwrap();
    //let mut sys = System::read("tests/diala.top");
    let top = Topology::read("tests/diala.top");
    let (mut pbc, mut coords) = trajectory::reader::read_coords("tests/diala.crd");
    let mut forces = vec![[0.0; 3]; coords.len()];
    let mut velocities = vec![[0.0; 3]; coords.len()];

    let atoms = top.get_atoms();
    let masses = atoms.iter().map(|a| a.mass).collect::<Vec<f32>>();
    let charges = atoms.iter().map(|a| a.charge).collect::<Vec<f32>>();

    let ffield = Forces::new(&top);
    let mut vv = VelocityVerlet::new(0.0000001, coords.len());
    let mut traj = TrajectoryWriter::new("tests/simple.traj", 1, false);
    let mut u = 0.0;
    (0..1000).for_each(|_| {
        u = vv.step(&ffield, &mut coords, &mut forces, &mut velocities, &masses);
        if vv.step % traj.interval == 0 {
            traj.write(&coords, &vv.time);
        }
    });
    println!("done!")
}
