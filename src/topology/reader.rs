use crate::ffield::*;

use super::Topology;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::rc::Rc;

#[derive(Debug)]
enum Section {
    MAIN,
    NONBOND,
    BOND,
    ANGLE,
    DIHEDRAL,
}

pub fn read(filename: &str) -> Topology {
    let mut top = Topology::new();
    let f = File::open(filename).unwrap();

    let mut section = Section::MAIN;
    let mut section_counter = 0;

    for line in BufReader::new(f).lines() {
        if let Ok(line) = line {
            if line.is_empty() {
                continue;
            }
            if line.starts_with("MAIN") {
                section = Section::MAIN;
                section_counter = 0;
            } else if line.starts_with("NONBOND") {
                section = Section::NONBOND;
                section_counter = 0;
            } else if line.starts_with("BOND") {
                section = Section::BOND;
                section_counter = 0;
            } else if line.starts_with("ANGLE") {
                section = Section::ANGLE;
                section_counter = 0;
            } else if line.starts_with("DIHEDRAL") {
                section = Section::DIHEDRAL;
                section_counter = 0;
            } else {
                section_counter += 1;
            }

            if section_counter > 0 {
                let fields: Vec<&str> = line.split_whitespace().collect();

                match section {
                    Section::MAIN => {
                        let atomtype = fields[0].to_string();
                        let name = fields[1].to_string();
                        let element = fields[2].to_string();
                        let mass = fields[3].parse::<f64>().unwrap();
                        let vdw = fields[4].parse::<f64>().unwrap();
                        let charge = fields[5].parse::<f64>().unwrap();
                        let x = fields[6].parse::<f64>().unwrap();
                        let y = fields[7].parse::<f64>().unwrap();
                        let z = fields[8].parse::<f64>().unwrap();
                        top.add_atom(atomtype, name, element, mass, vdw, charge, [x, y, z]);
                    }
                    Section::NONBOND => {
                        let a = fields[0].parse::<usize>().unwrap();
                        let c12 = fields[1].parse::<f64>().unwrap();
                        let c6 = fields[2].parse::<f64>().unwrap();
                        top.atoms[a].borrow_mut().LJ.c12 = c12;
                        top.atoms[a].borrow_mut().LJ.c6 = c6;

                    }
                    Section::BOND => {
                        let a1 = fields[0].parse::<usize>().unwrap();
                        let a2 = fields[1].parse::<usize>().unwrap();
                        let k = fields[2].parse::<f64>().unwrap();
                        let r0 = fields[3].parse::<f64>().unwrap();
                        top.add_bond(a1, a2);
                        top.ff.bond.push(Bond::new(
                            k,
                            r0,
                            [Rc::clone(&top.atoms[a1]), Rc::clone(&top.atoms[a2])],
                        ));
                    }
                    Section::ANGLE => {
                        let a1 = fields[0].parse::<usize>().unwrap();
                        let a2 = fields[1].parse::<usize>().unwrap();
                        let a3 = fields[2].parse::<usize>().unwrap();
                        let k = fields[3].parse::<f64>().unwrap();
                        let t0 = fields[4].parse::<f64>().unwrap();
                        top.ff.angle.push(Angle::new (
                            k,
                            t0,
                            [
                                Rc::clone(&top.atoms[a1]),
                                Rc::clone(&top.atoms[a2]),
                                Rc::clone(&top.atoms[a3]),
                            ],
                       ));
                    }
                    Section::DIHEDRAL => {
                        let a1 = fields[0].parse::<usize>().unwrap();
                        let a2 = fields[1].parse::<usize>().unwrap();
                        let a3 = fields[2].parse::<usize>().unwrap();
                        let a4 = fields[3].parse::<usize>().unwrap();
                        let k = fields[4].parse::<f64>().unwrap();
                        let t0 = fields[5].parse::<f64>().unwrap();
                        top.ff.dihedral.push(Dihedral::new (
                            k,
                            t0,
                            [
                                Rc::clone(&top.atoms[a1]),
                                Rc::clone(&top.atoms[a2]),
                                Rc::clone(&top.atoms[a3]),
                                Rc::clone(&top.atoms[a4]),
                            ],
                        ));
                    }
                }
            }
        }
    }
    return top;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_reads() {
        let top = read("tests/simple.topol");
        top.atoms[0].borrow_mut().name = "AAAAA".to_string();

        println!("{:#?}", top)
    }
}
