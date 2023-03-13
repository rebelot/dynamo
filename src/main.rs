use dynamo::integrator::verlet::VelocityVerlet;
use dynamo::system::System;


fn main() {
    let sys = &mut System::read("tests/simple.sys");
    let vv = &mut VelocityVerlet::new(0.001);
    while vv.time < 1_000.0 {
        vv.step(sys);
    }
    println!("{:#?}", *sys)
}
