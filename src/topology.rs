pub mod atom;
pub mod molecule;
mod reader;
use atom::*;
use molecule::*;
use std::collections::HashMap;

use crate::ffield::Forces;

#[derive(Debug, Default)]
pub struct Defaults {
    nb_func: String,
    comb_rule: String,
    ljscale: Option<f32>,
    qqscale: Option<f32>,
}

#[derive(Debug, Default)]
pub struct Topology {
    pub atomtypes: HashMap<String, AtomTypeParams>,
    pub molecules: Vec<Molecule>,
    pub natoms: usize,
    pub nmols: usize,
    pub pbc: [f32; 3],
    pub defaults: Defaults,
}

impl Topology {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_defaults(
        &mut self,
        nb_func: String,
        comb_rule: String,
        ljscale: Option<f32>,
        qqscale: Option<f32>,
    ) {
        self.defaults.nb_func = nb_func;
        self.defaults.comb_rule = comb_rule;
        self.defaults.ljscale = ljscale;
        self.defaults.qqscale = qqscale;
    }

    pub fn add_atomtype(&mut self, atomtype: String, element: u32, mass: f32, v: f32, w: f32) {
        let params = AtomTypeParams {
            element,
            mass,
            v,
            w,
        };
        self.atomtypes.insert(atomtype, params);
    }

    pub fn add_bonded_interaction(&mut self, funct: &str, params: Vec<f32>) {
        self.molecules[self.nmols - 1]
            .bonded_interactions
            .push((funct.to_owned(), params))
    }

    pub fn add_molecule(&mut self, name: String, nmols: usize, nbexc: usize) {
        let molecule = Molecule::new(self.nmols, name, nmols, nbexc);
        self.molecules.push(molecule);
        self.nmols += 1;
    }

    pub fn add_exclpairs(&mut self, i: usize, excl: Vec<usize>) {
        let mut excl = excl;
        excl.sort();
        self.molecules[self.nmols - 1].atoms[i]
            .excluded
            .extend_from_slice(&excl);
    }

    pub fn add_atom(
        &mut self,
        atomtype: String,
        name: String,
        resnum: usize,
        resname: String,
        charge: f32,
    ) {
        if !self.atomtypes.contains_key(&atomtype) {
            panic!("Undefined atom type: {}", atomtype);
        }
        let atom = atom::Atom::new(
            self.natoms,
            atomtype.to_owned(),
            name,
            resnum,
            resname,
            self.atomtypes[&atomtype].element,
            self.atomtypes[&atomtype].mass,
            charge,
            self.atomtypes[&atomtype].v,
            self.atomtypes[&atomtype].w,
        );
        self.molecules[self.nmols - 1].atoms.push(atom);
        self.natoms += 1;
    }

    pub fn read(filename: &str) -> Self {
        reader::read(filename)
    }

    pub fn build_forces(&mut self) -> Forces {
        let mut forces = Forces::new();
        self.molecules.iter().for_each(|mol| {
            (0..mol.nmols).for_each(|i| {
                mol.bonded_interactions.iter().for_each(|interaction| {
                    forces.parse_interaction(interaction.to_owned(), mol.atoms.len() * i)
                })
            })
        });
        forces
    }
}
