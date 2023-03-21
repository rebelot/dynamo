use dynamo::integrator::verlet::VelocityVerlet;
use dynamo::system::System;
use dynamo::trajectory::writer::TrajectoryWriter;

fn main() {
    let mut sys = System::read("tests/simple.sys");
    let mut vv = VelocityVerlet::new(0.001);
    let mut traj = TrajectoryWriter::new("tests/simple.traj", 1, false);
    while vv.time < 1_000.0 {
        vv.step(&mut sys);
        traj.write(&sys, &vv.step, &vv.time);
    }
    // println!("{:#?}", sys)
}
