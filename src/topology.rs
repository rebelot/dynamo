pub mod atom;
pub mod molecule;
mod reader;
use atom::*;
use molecule::*;
use std::collections::HashMap;

use crate::ffield::Forces;

#[derive(Debug, Default)]
pub struct Defaults {
    pub nb_func: String,
    pub comb_rule: String,
    pub ljscale: Option<f32>,
    pub qqscale: Option<f32>,
}

#[derive(Debug, Default)]
pub struct Topology {
    pub atomtypes: HashMap<String, AtomTypeParams>,
    pub molecules: Vec<Molecule>,
    natoms: usize,
    nmols: usize,
    pub defaults: Defaults,
}

impl Topology {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_defaults(
        &mut self,
        nb_func: &str,
        comb_rule: &str,
        ljscale: Option<f32>,
        qqscale: Option<f32>,
    ) {
        self.defaults.nb_func = nb_func.to_string();
        self.defaults.comb_rule = comb_rule.to_string();
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

    pub fn add_bonded_interaction(&mut self, moli: usize, interaction: &str) {
        self.molecules[moli].add_bonded_interaction(interaction);
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
        moli: usize,
        atomtype: &str,
        name: &str,
        resnum: usize,
        resname: &str,
        charge: f32,
    ) {
        if !self.atomtypes.contains_key(atomtype) {
            panic!("Undefined atom type: {}", atomtype);
        }
        let atom = atom::Atom::new(
            self.natoms,
            atomtype.to_owned(),
            name.to_owned(),
            resnum,
            resname.to_owned(),
            self.atomtypes[atomtype].element,
            self.atomtypes[atomtype].mass,
            charge,
            self.atomtypes[atomtype].v,
            self.atomtypes[atomtype].w,
        );
        self.molecules[moli].atoms.push(atom);
        self.natoms += 1;
    }

    pub fn get_atoms(&self) -> Vec<Atom> {
        let mut atoms = Vec::new();
        for mol in &self.molecules {
            for _ in 0..mol.nmols {
                for atom in &mol.atoms {
                    atoms.push(atom.clone());
                }
            }
        }
        atoms
    }

    pub fn read(filename: &str) -> Self {
        reader::parse(filename)
    }
}
