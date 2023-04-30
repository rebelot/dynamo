use dynamo::integrator::verlet::VelocityVerlet;
use dynamo::system::System;
use dynamo::trajectory::writer::TrajectoryWriter;

fn main() {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global().unwrap();
    // let mut sys = System::read("tests/simple.top");
    let mut vv = VelocityVerlet::new(0.001);
    let mut traj = TrajectoryWriter::new("tests/simple.traj", 1, false);
    (0..10_000).for_each(|_| {
        // vv.step(&mut sys);
        // traj.write(&sys, &vv.step, &vv.time);
    });
}
