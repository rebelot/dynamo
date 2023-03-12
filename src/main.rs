use std::{cell::RefCell, rc::Rc};

use dynamo::topology::{self, atom::Atom};

fn main() {
    let a = Atom::new(
        0,
        "x".to_owned(),
        "x".to_owned(),
        "x".to_owned(),
        0.0,
        0.0,
        0.0,
        [0.0, 0.0, 0.0],
    );
    add_one(&mut a.borrow_mut());
    println!("a: {:?}", a);
}

fn add_one(a: &mut Atom) {
    a.index += 1;
}
