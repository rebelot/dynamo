use super::Topology;
use std::fs::File;
use std::io::{BufRead, BufReader};

// Specification of the topology file format:
//  * Section MAIN
//      nb_func comb_rule [ljscale qqscale]
//  * Section ATOMTYPES
//      type element mass charge c0 c1
//  * Section MOL
//      name nbexc
//  * Section ATOMS
//      index type name resnum resname [charge c0 c1]
//  * Section BONDED
//      funct ...params
//  * Section
//      a b c

enum Section {
    Main,
    AtomTypes,
    Molecule,
    Atoms,
    Bonds,
    ExclPairs,
    System,
}

impl Default for Section {
    fn default() -> Self {
        Self::Main
    }
}

pub fn parse(filename: &str) -> Topology {
    let mut top = Topology::new();

    let f = File::open(filename).unwrap();

    let mut section: Section = Default::default();
    let mut section_counter = 0;

    for line in BufReader::new(f).lines().map_while(Result::ok) {
        if line.is_empty() {
            continue;
        }
        if line.starts_with("MAIN") {
            section = Section::Main;
            section_counter = 0;
        } else if line.starts_with("ATOMTYPES") {
            section = Section::AtomTypes;
            section_counter = 0;
        } else if line.starts_with("MOL") {
            section = Section::Molecule;
            section_counter = 0;
        } else if line.starts_with("ATOMS") {
            section = Section::Atoms;
            section_counter = 0;
        } else if line.starts_with("BONDS") {
            section = Section::Bonds;
            section_counter = 0;
        } else if line.starts_with("EXCLPAIRS") {
            section = Section::ExclPairs;
            section_counter = 0;
        } else if line.starts_with("SYSTEM") {
            section = Section::System;
            section_counter = 0;
        } else {
            section_counter += 1;
        }

        if section_counter > 0 {
            let fields: Vec<&str> = line.split_whitespace().collect();

            match section {
                Section::Main => top.set_defaults(
                    fields[0],
                    fields[1],
                    fields[2].parse::<f32>().ok(),
                    fields[2].parse::<f32>().ok(),
                ),

                Section::AtomTypes => top.add_atomtype(
                    fields[0].to_owned(),
                    fields[1].parse::<u32>().unwrap(),
                    fields[2].parse::<f32>().unwrap(),
                    fields[3].parse::<f32>().unwrap(),
                    fields[4].parse::<f32>().unwrap(),
                ),

                Section::Molecule => {
                    top.add_molecule(fields[0].to_owned(), 0, fields[1].parse::<usize>().unwrap())
                }

                Section::Atoms => top.add_atom(
                    top.nmols - 1,
                    fields[1],
                    fields[2],
                    fields[3].parse::<usize>().unwrap(),
                    fields[4],
                    fields[5].parse::<f32>().unwrap(),
                ),

                Section::ExclPairs => top.add_exclpairs(
                    fields[0].parse::<usize>().unwrap(),
                    fields[1..]
                        .iter()
                        .map(|a| a.parse::<usize>().unwrap())
                        .collect(),
                ),

                Section::System => {
                    let molname = fields[0].parse::<String>().unwrap();
                    for mol in top.molecules.iter_mut() {
                        if mol.name == molname {
                            mol.nmols = fields[1].parse::<usize>().unwrap();
                            break;
                        }
                    }
                }

                Section::Bonds => top.add_bonded_interaction(top.nmols - 1, &fields.join(" ")),
            }
        }
    }
    top
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_reads() {
        let sys = parse("tests/diala.top");
        println!("{:#?}", sys);
    }
}
