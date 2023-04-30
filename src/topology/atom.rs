#[derive(Debug)]
pub struct AtomTypeParams {
    pub element: u32,
    pub mass: f32,
    pub v: f32,
    pub w: f32,
}

#[derive(Debug, Clone)]
// Stores particle information
pub struct Atom {
    pub index: usize,
    pub resname: String,
    pub resnum: usize,
    pub atomtype: String,
    pub name: String,
    pub element: u32,
    pub mass: f32,
    pub charge: f32,
    pub v: f32,
    pub w: f32,
    pub excluded: Vec<usize>,
}

impl Atom {
    pub fn new(
        index: usize,
        atomtype: String,
        name: String,
        resnum: usize,
        resname: String,
        element: u32,
        mass: f32,
        charge: f32,
        v: f32,
        w: f32,
    ) -> Atom {
        Atom {
            index,
            resnum,
            resname,
            atomtype,
            name,
            element,
            mass,
            charge,
            v,
            w,
            excluded: Vec::new(),
        }
    }
}
