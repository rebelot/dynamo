use super::atom::*;

#[derive(Debug, Clone)]
pub struct Molecule {
    pub index: usize,
    pub name: String,
    pub nmols: usize,
    pub nbexc: usize,
    pub bonded_interactions: Vec<(String, Vec<f32>)>,
    pub atoms: Vec<Atom>,
}

impl Molecule {
    pub fn new(index: usize, name: String, nmols: usize, nbexc: usize) -> Molecule {
        Molecule {
            index,
            name,
            nmols,
            nbexc,
            atoms: Vec::new(),
            bonded_interactions: Vec::new()
        }
    }
}
