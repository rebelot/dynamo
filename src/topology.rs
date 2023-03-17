pub mod atom;

#[derive(Debug)]
pub struct Topology {
    pub atoms: Vec<atom::Atom>,
}

impl Topology {
    pub fn new() -> Self {
        return Topology { atoms: Vec::new() };
    }

    pub fn add_atom(
        &mut self,
        atomtype: String,
        name: String,
        element: String,
        mass: f64,
        vdw: f64,
        charge: f64,
        coords: [f64; 3],
    ) -> usize {
        let index = self.atoms.len();
        let atom = atom::Atom::new(index, atomtype, name, element, mass, vdw, charge, coords);
        self.atoms.push(atom);
        return index;
    }

    pub fn init_forces(&mut self) {
        for a in self.atoms.iter_mut() {
            for i in 0..3 {
                a.prev_force[i] = a.force[i];
                a.force[i] = 0.0;
            }
        }
    }
}
